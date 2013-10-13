#include "Python.h"
#include "svm.h"

#define MAX_LINE_LENGTH 2048

struct svm_model* Model;
double DecisionValues[10]; // hacky
struct svm_node* SVMFeatures;
int SVMFeatureAllocation;

// Assume no more than 128 features!
double MinFeatureValue[128];
double MaxFeatureValue[128];

static PyObject* ex_foo(PyObject* self, PyObject* args)
{
    printf("Hello, world\n");
    Py_INCREF(Py_None);
    return Py_None;
}

// Copy one line (up to a \r or \n character) from a source buffer to a target buffer.
// Optionally, strip out spaces.  Return the position just AFTER the end of the line.
// (If a line ends in \r\n, we'll end up processing the line, and then one empty line; that's okay)
// If a line is very long, we stop copying, and skip over the rest of it.
int CopyBufferLine(char* Source, int BufferPos, int BufferEnd, char* LineBuffer, int StripWhitespace)
{
    int LinePos = 0;
    int LineComplete = 0;
    int Chars = 0;
    int Skipping = 0;
    //
    while (!LineComplete)
    {
        if (BufferPos > BufferEnd)
        {
            // Our line extends off the edge of the buffer.  That's probably a Bad Thing.
            printf("** Warning: Ran off the edge of the buffer in CopyBufferLine.  Line too ling?\n");
            LineBuffer[LinePos] = '\0';
            return BufferPos;
        }
        switch (Source[BufferPos])
        {
        case ' ':
            if (StripWhitespace)
            {
                BufferPos++;
            }
            else
            {
                if (!Skipping)
                {
                    LineBuffer[LinePos++] = Source[BufferPos];
                }
                BufferPos++;
                Chars++;
            }
            break;
        case '\r':
        case '\n':
            LineBuffer[LinePos] = '\0';
            BufferPos++;
            LineComplete = 1;
            break;
        case '\0':
            LineBuffer[LinePos] = '\0';
            LineComplete = 1;
            break;
        default:
            if (!Skipping)
            {
                LineBuffer[LinePos++] = Source[BufferPos];
            }
            BufferPos++;
            Chars++;
            break;
        }
        if (Chars == MAX_LINE_LENGTH - 1)
        {
            printf("** Error: Line too long!  Truncating line.");
            // Read the rest of the chars, but don't write them:
            Chars = 0;
            Skipping = 1;
        }
    }
    return BufferPos;
}

static PyObject* PyLoadScaling(PyObject* self, PyObject* args)
{
    char* FilePath;
    char* FileText;
    char LineBuffer[MAX_LINE_LENGTH];
    FILE* ScaleFile;
    int BufferPos; 
    int BufferEnd;
    int FeatureNumber;
    char* StrValue;
    double MinValue;
    double MaxValue;
    //
    if (!PyArg_ParseTuple(args, "s", &FilePath))
    {
        return (PyObject*)-1; // Return -1 to signal that the object can't be created
    }
    ScaleFile = fopen(FilePath, "rb");
    if (!ScaleFile)
    {
        printf("** Error: Can't open file '%s'\n", FilePath);
        return (PyObject*)-1; // Return -1 to signal that the object can't be created
    }
    FileText = (char*)calloc(sizeof(char), 10240);
    BufferEnd = fread(FileText, sizeof(char), 10240, ScaleFile);
    BufferPos = 0;
    while (1)
    {
        if (!FileText[BufferPos])
        {
            break;
        }
        BufferPos = CopyBufferLine(FileText, BufferPos, BufferEnd, LineBuffer, 0);
        //printf("Line parsed: '%s'\n", LineBuffer);
        StrValue = strtok(LineBuffer, " \t");
        if (!StrValue)
        {
            continue;
        }
        FeatureNumber = atoi(StrValue);
        if (FeatureNumber <= 0)
        {
            continue;
        }
        StrValue = strtok(NULL, " \t");
        if (!StrValue)
        {
            continue;
        }
        MinValue = atof(StrValue);
        StrValue = strtok(NULL, " \t");
        if (!StrValue)
        {
            continue;
        }
        MaxValue = atof(StrValue);
        MinFeatureValue[FeatureNumber - 1] = MinValue;
        MaxFeatureValue[FeatureNumber - 1] = MaxValue;
        //printf("Feature %d: Range %f...%f\n", FeatureNumber, MinValue, MaxValue);
    }
    fclose(ScaleFile);
    Py_INCREF(Py_None);
    return Py_None;
}


void ScaleSVMFeatures(struct svm_node* Features, int FeatureCount)
{
    int FeatureIndex;
    double Value;
    double Range;
    //
    for (FeatureIndex = 0; FeatureIndex < FeatureCount; FeatureIndex++)
    {
        Value = Features[FeatureIndex].value;
        Range = MaxFeatureValue[FeatureIndex] - MinFeatureValue[FeatureIndex];
        if (Value <= MinFeatureValue[FeatureIndex])
        {
            SVMFeatures[FeatureIndex].value = -1.0;
            continue;
        }
        if (Value >= MaxFeatureValue[FeatureIndex])
        {
            Features[FeatureIndex].value = 1.0;
            continue;
        }
        Features[FeatureIndex].value = -1.0 + 2.0 * (Value - MinFeatureValue[FeatureIndex]) / Range;
    }
}

static PyObject* PyLoadModel(PyObject* self, PyObject* args)
{
    char* FilePath;
    if (!PyArg_ParseTuple(args, "s", &FilePath))
    {
        return (PyObject*)-1; // Return -1 to signal that the object can't be created
    }

    // Free the old model, if any:
    if (Model)
    {
        svm_destroy_model(Model);
        Model = NULL;
    }

    // Load model from specified file:
    Model = svm_load_model(FilePath);

    Py_INCREF(Py_None);
    return Py_None;
}

PyObject* PyScoreHelper(PyObject* FeatureList, int ScaleFlag)
{
    int FeatureIndex;
    int FeatureCount;
    int SequenceType;
    //
    if (PyList_Check(FeatureList))
    {
        SequenceType = 1;
        FeatureCount = PyList_Size(FeatureList);
    }
    else if PyTuple_Check(FeatureList)
    {
        SequenceType = 2;
        FeatureCount = PyTuple_Size(FeatureList);
    }
    else
    {
        printf("** Error in PyScore vector: Illegal argument (not a vector or tuple)\n");
        return (PyObject*)-1; // Return -1 to signal that the object can't be created
    }

    // Allocate SVMFeatures, if necessary:
    if (FeatureCount >= SVMFeatureAllocation)
    {
        if (SVMFeatures)
        {
            free(SVMFeatures);
        }
        SVMFeatures = (struct svm_node*)malloc((FeatureCount + 1)  * sizeof(struct svm_node));
        SVMFeatureAllocation = FeatureCount + 1;
    }

    // Populate SVMFeatures:
    for (FeatureIndex = 0; FeatureIndex < FeatureCount; FeatureIndex++)
    {
        SVMFeatures[FeatureIndex].index = FeatureIndex + 1;
        if (SequenceType == 1)
        {
            SVMFeatures[FeatureIndex].value = PyFloat_AsDouble(PyList_GetItem(FeatureList, FeatureIndex));
        }
        else
        {
            SVMFeatures[FeatureIndex].value = PyFloat_AsDouble(PyTuple_GetItem(FeatureList, FeatureIndex));
        }
    }
    if (ScaleFlag)
    {
        ScaleSVMFeatures(SVMFeatures, FeatureCount);
    }
    SVMFeatures[FeatureCount].index = -1;
    // Predict, and return:
    svm_predict_values(Model, SVMFeatures, DecisionValues);
    return PyFloat_FromDouble(DecisionValues[0]);
}

static PyObject* PyScoreVector(PyObject* self, PyObject* args)
{
    PyObject* FeatureList;
    //
    if (!PyArg_ParseTuple(args, "O", &FeatureList))
    {
        return (PyObject*)-1; // Return -1 to signal that the object can't be created
    }
    return PyScoreHelper(FeatureList, 0);
}

static PyObject* PyScaleAndScoreVector(PyObject* self, PyObject* args)
{
    PyObject* FeatureList;
    //
    if (!PyArg_ParseTuple(args, "O", &FeatureList))
    {
        return (PyObject*)-1; // Return -1 to signal that the object can't be created
    }
    return PyScoreHelper(FeatureList, 1);
}

static PyMethodDef PySVM_methods[] = {
    {"foo", ex_foo, METH_VARARGS, "foo() doc string"},
    {"LoadModel", PyLoadModel, METH_VARARGS, "Load an SVM model from disk"},
    {"Score", PyScoreVector, METH_VARARGS, "Score a (pre-scaled) vector"},
    {"LoadScaling", PyLoadScaling, METH_VARARGS, "Load feature scaling parameters from a file"},
    {"ScaleAndScore", PyScaleAndScoreVector, METH_VARARGS, "Scale and score a feature-vector"},
    {NULL, NULL}
};

PyMODINIT_FUNC initPySVM(void)
{
    Py_InitModule("PySVM", PySVM_methods);
    Model = NULL;
    SVMFeatures = NULL;
    SVMFeatureAllocation = 0;
}
