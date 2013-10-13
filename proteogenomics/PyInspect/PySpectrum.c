// PySpectrum: Python wrapper for an MSSpectrum object.
#include "CMemLeak.h"
#include "PySpectrum.h"
#include "PyUtils.h"
#include "Trie.h"
#include "Tagger.h"
#include "Score.h"
#include "Mods.h"
#include "CMemLeak.h"
#include "FreeMod.h"
#include "BN.h"
#include "ChargeState.h"
#include "Scorpion.h"
#include "SVM.h"
#include "IonScoring.h"
#include "Errors.h"

// Important note: These type objects MUST be defined here in the .c file, not in the header! 
// Otherwise, a copy is built for each including file, and these copies are not 
// updated during module initialization.  (The result is that MSPeak objects instantiated
// from code don't have their members set right, so their attributes can't be accessed)
PyTypeObject PySpectrumType = 
{
    PyObject_HEAD_INIT(NULL)
    0, //ob_size
    "PyInspect.PySpectrum",  //tp_name
    sizeof(PySpectrum),  //tp_basicsize
    0,                         //tp_itemsize
    PySpectrumDealloc,                         //tp_dealloc
    0,                         //tp_print
    0,                         //tp_getattr
    0,                         //tp_setattr
    0,                         //tp_compare
    0,                         //tp_repr
    0,                         //tp_as_number
    0,                         //tp_as_sequence
    0,                         //tp_as_mapping
    0,                         //tp_hash 
    0,                         //tp_call
    0,                         //tp_str
    0,                         //tp_getattro
    0,                         //tp_setattro
    0,                         //tp_as_buffer
    Py_TPFLAGS_DEFAULT,        //tp_flags
    "MS spectrum",           // tp_doc 
    0,                       // tp_traverse 
    0,                       // tp_clear 
    0,                       // tp_richcompare 
    0,                       // tp_weaklistoffset 
    0,                       // tp_iter 
    0,                       // tp_iternext 
    PySpectrumMethods,             // tp_methods 
    PySpectrumMembers,             // tp_members 
    PySpectrumGetSet,           // tp_getset 
    0,                         // tp_base 
    0,                         // tp_dict 
    0,                         // tp_descr_get 
    0,                         // tp_descr_set 
    0,                         // tp_dictoffset 
    (initproc)PySpectrumInit,      // tp_init 
    0,                         // tp_alloc 
    PySpectrumNew,                 // tp_new 
};

extern PRMBayesianModel* InteractiveBN; // lives in PyInspect.c

TrieTag* TagGraphGenerateTags(TagGraph* Graph, MSSpectrum* Spectrum, int* TagCount, 
    int MaximumTagCount, SpectrumTweak* Tweak, float TagEdgeScoreMultiplier,
    PRMBayesianModel* Model);

// __new__ method of PySpectrum; call this in C code to create new PySpectrum objects.
// (It's expected that PySpectrumInit gets called too)
PyObject* PySpectrumNew(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PySpectrum* self;
    //
    self = (PySpectrum*)type->tp_alloc(type, 0);
    if (self != NULL) 
    {
        // Perform non-parameterized initialization here. 
        //memset(self->NumberToMSIndex, -1, MAX_MS_SCAN);
    }
    return (PyObject*)self;
}

// Called when Inspect.Spectrum() is instantiated from Python code.
// Parse the specified spectrum file!
PyObject* PySpectrumInit(PySpectrum* self, PyObject* args, PyObject* kwds)
{
    char* FilePath;
    int FilePosition = 0; // Default: byte offset 0
    SpectrumNode* Node;
    FILE* SpectrumFile;
    int LoadResult;
    //
    // Constructor argument: the path to a mass-spec run output file
    if (!PyArg_ParseTuple(args, "s|i", &FilePath, &FilePosition))
    {
        return (PyObject*)-1;
    }
    Node = (SpectrumNode*)calloc(1, sizeof(SpectrumNode));
    Node->FilePosition = FilePosition;
    Node->ScanNumber = 0;
    Node->InputFile = (InputFileNode*)calloc(1, sizeof(InputFileNode));
    strncpy(Node->InputFile->FileName, FilePath, MAX_FILENAME_LEN);
    strncpy(self->FileName, FilePath, MAX_FILENAME_LEN);
    // Guess the file format:
    Node->InputFile->Format = GuessSpectrumFormatFromExtension(FilePath);

    Node->Spectrum = (MSSpectrum*)calloc(1, sizeof(MSSpectrum));
    Node->Spectrum->Node = Node;

    // Special case: If it's a .ms2 extension, it could be "colon format" or standard MS2:
    if (Node->InputFile->Format == SPECTRUM_FORMAT_MS2_COLONS)
    {
        Node->InputFile->Format = GuessSpectrumFormatFromHeader(FilePath, Node->Spectrum);
    }
    SpectrumFile = fopen(FilePath, "rb");
    if (!SpectrumFile)
    {
        sprintf(PythonErrorString, "** Error: Unable to open spectrum file '%s'\n", FilePath);
        ReportPythonError();
        // In an Init function, we must return -1 to indicate that an object can't
        // be created.  (Normally we return NULL for failure!)
        return (PyObject*)-1;
    }
    fseek(SpectrumFile, Node->FilePosition, 0);
    LoadResult = SpectrumLoadFromFile(Node->Spectrum, SpectrumFile);
    fclose(SpectrumFile);
    if (!LoadResult)
    {
        sprintf(PythonErrorString, "** Error: Unable to parse spectrum from %s:%d\n", FilePath, Node->FilePosition);
        ReportPythonError();
        // In an Init function, we must return -1 to indicate that an object can't
        // be created.  (Normally we return NULL for failure!)
        return (PyObject*)-1;
    }
    WindowFilterPeaks(Node->Spectrum, 0, 0);
    IntensityRankPeaks(Node->Spectrum);
    self->Spectrum = Node->Spectrum;
    TweakSpectrum(Node);
    return 0;
}

int GuessCharge(SpectrumNode* Node, int MatchMass)
{
    int BestDiff = -1;
    int BestCharge;
    int Charge;
    int Mass;
    int Diff;
    //
    for (Charge = 1; Charge < 5; Charge++)
    {
        Mass = Node->Spectrum->MZ * Charge - (Charge - 1) * HYDROGEN_MASS;
        Diff = abs(MatchMass - Mass);
        if (BestDiff < 0 || Diff < BestDiff)
        {
            BestDiff = Diff;
            BestCharge = Charge;
        }
    }
    return BestCharge;
}

// Helper function for both LabelPeaks() and Score(); the code paths overlap
PyObject* ScoreHelper(PySpectrum* self, char* PeptideString, int Charge, int VerboseFlag, int LabelPeaksFlag, int ReturnScoreDetails)
{
    Peptide* Match;
    float MQScore;
    SpectrumNode* Node;
    int PeakIndex;
    SpectralPeak* Peak;
    PyObject* LabeledPeakList;
    PyObject* LabeledPeakTuple;
    PyListObject* List;
    int FeatureIndex;
    //

    Node = self->Spectrum->Node;
    Match = GetPeptideFromAnnotation(PeptideString);
    if (!Match)
    {
        sprintf(PythonErrorString, "** Error: Unable to parse peptide annotation '%s'\n", PeptideString);
        ReportPythonError();
        Py_INCREF(Py_None);
        return Py_None;
    }
    if (!Charge)
    {
        Charge = GuessCharge(Node, Match->ParentMass);    
    }
    Node->Tweaks[0].ParentMass = Match->ParentMass;
    Node->Tweaks[0].Charge = Charge;
    Match->Tweak = Node->Tweaks;
    Node->Spectrum->Charge = Charge;
    GlobalOptions->DigestType = DIGEST_TYPE_TRYPSIN;
    ComputeMQScoreFeatures(Node->Spectrum, Match, Match->ScoreFeatures, VerboseFlag);

#ifdef MQSCORE_USE_SVM
    MQScore = SVMComputeMQScore(Node->Spectrum, Match, Match->ScoreFeatures);
#else
    MQScore = LDAComputeMQScore(Node->Spectrum, Match, Match->ScoreFeatures);
#endif

    if (VerboseFlag)
    {
        // Print out the ion-type categorization of all peaks:
        printf("\n");
        printf("Score %s on spectrum %s:%d\n", PeptideString, Node->InputFile->FileName, Node->FilePosition);
        for (PeakIndex = 0; PeakIndex < Node->Spectrum->PeakCount; PeakIndex++)
        {
            Peak = Node->Spectrum->Peaks + PeakIndex;
            printf("%.2f\t%.2f\t%s\t%d\t\n", Peak->Mass / (float)MASS_SCALE, Peak->Intensity, GetFragmentTypeName(Peak->IonType), Peak->AminoIndex);
        }
    }
    
    if (LabelPeaksFlag)
    {
        // Return a list of peaks.
        LabeledPeakList = PyList_New(0);
        for (PeakIndex = 0; PeakIndex < Node->Spectrum->PeakCount; PeakIndex++)
        {
            Peak = Node->Spectrum->Peaks + PeakIndex;
            LabeledPeakTuple = Py_BuildValue("ffsi", Peak->Mass / (float)MASS_SCALE, Peak->Intensity, GetFragmentTypeName(Peak->IonType), Peak->AminoIndex);
            PyList_Append(LabeledPeakList, LabeledPeakTuple);
            // The call to PyList_Append has incremented the refcount of the 
            // tuple to 2.  We're abandoning our reference to the tuple now, 
            // so we decrease its reference count:
            Py_DECREF(LabeledPeakTuple);
            LabeledPeakTuple = NULL; // just to be explicit about it!
        }
        FreePeptideNode(Match);
        return LabeledPeakList;
    }
    else
    {
        if (ReturnScoreDetails)
        {
            List = (PyListObject*)PyList_New(0);
            PyList_Append((PyObject*)List, PyFloat_FromDouble(MQScore));
            for (FeatureIndex = 0; FeatureIndex < 7; FeatureIndex++)
            {
                PyList_Append((PyObject*)List, PyFloat_FromDouble(Match->ScoreFeatures[FeatureIndex]));
            }
            FreePeptideNode(Match);
            //Py_DECREF(List); // Important to do this!
            return (PyObject*)List;
        }
        else
        {
            FreePeptideNode(Match);
            return PyFloat_FromDouble(MQScore);
        }
    }
}

PyObject* PySpectrumLabelPeaks(PySpectrum* self, PyObject* args, PyObject* kwargs)
{
    char* PeptideString;
    int Charge = 0;
    int VerboseFlag = 0;
    // 
    if (!PyArg_ParseTuple(args, "s|ii", &PeptideString, &Charge, &VerboseFlag))
    {
        return NULL;
    }
    return ScoreHelper(self, PeptideString, Charge, VerboseFlag, 1, 0);
}

PyObject* PySpectrumScore(PySpectrum* self, PyObject* args, PyObject* kwargs)
{
    char* PeptideString;
    int Charge = 0;
    int VerboseFlag = 0;
    if (!PyArg_ParseTuple(args, "s|ii", &PeptideString, &Charge, &VerboseFlag))
    {
        return NULL;
    }
    return ScoreHelper(self, PeptideString, Charge, VerboseFlag, 0, 0);
}

PyObject* PySpectrumScoreDetailed(PySpectrum* self, PyObject* args, PyObject* kwargs)
{
    char* PeptideString;
    int Charge = 0;
    int VerboseFlag = 0;
    if (!PyArg_ParseTuple(args, "s|ii", &PeptideString, &Charge, &VerboseFlag))
    {
        return NULL;
    }
    return ScoreHelper(self, PeptideString, Charge, VerboseFlag, 0, 1);
}


// Deallocate an PySpectrum object.  
void PySpectrumDealloc(PyObject* selfobject)
{
    PySpectrum* self = (PySpectrum*)selfobject;
    if (self)
    {
        if (self->Spectrum)
        {
            // The PySpectrum object wraps a Spectrum object, but also a SpectrumNode and an InputFileNode.  
            // So, free those as well:
            if (self->Spectrum->Node->InputFile)
            {
                free(self->Spectrum->Node->InputFile);
                self->Spectrum->Node->InputFile = NULL;
            }
            if (self->Spectrum->Node)
            {
                FreeSpectrumNode(self->Spectrum->Node);
            }
            else
            {
                FreeSpectrum(self->Spectrum);
            }
            if (self->MatchList)
            {
                Py_DECREF(self->MatchList);
            }
        }
        self->ob_type->tp_free((PyObject*)self);
    }
}

PyObject* PySpectrumGetPeakCount(PySpectrum* self, PyObject* args, PyObject* kwargs)
{
    if (!PyArg_ParseTuple(args, ""))
    {
        return NULL;
    }
    return Py_BuildValue("i", self->Spectrum->PeakCount);
}

//PyObject* PySpectrumGetCharge(PySpectrum* self, PyObject* args, PyObject* kwargs)
//{
//    if (!PyArg_ParseTuple(args, ""))
//    {
//        return NULL;
//    }
//    return Py_BuildValue("i", self->Spectrum->Charge);
//}


// Run parent mass correction.  Return a list of tuples of the 
// form (Mass, ModelScore).
PyObject* PySpectrumCorrectParentMass(PySpectrum* self, PyObject* args, PyObject* kwargs)
{
    int CorrectChargeTemp = -1;
    int PMToleranceTemp = -1;
    int VerboseFlag = 0;
    PMCSpectrumInfo* SpectrumInfo;
    PyObject* FeatureTupleList;
    PyObject* FeatureTuple;
    PMCInfo* Info;
    MSSpectrum* Spectrum = self->Spectrum;

    //
    if (!PyArg_ParseTuple(args, "|ii", &PMToleranceTemp, &CorrectChargeTemp))
    {
        return NULL;
    }

    SpectrumInfo = GetPMCSpectrumInfo(Spectrum);
    PerformPMC(SpectrumInfo);
    
    FeatureTupleList = PyList_New(0);
    for (Info = SpectrumInfo->Head; Info; Info = Info->Next)
    {
        FeatureTuple = PyTuple_New(2);
        PyTuple_SetItem(FeatureTuple, 0, PyFloat_FromDouble(Info->ParentMass / (float)MASS_SCALE));
        PyTuple_SetItem(FeatureTuple, 1, PyFloat_FromDouble(Info->SVMScore));
        PyList_Append(FeatureTupleList, FeatureTuple);
    }

    //return ScoreHelper(self, PeptideString, Charge, VerboseFlag, 0, 0);
    FreePMCSpectrumInfo(SpectrumInfo);
    return FeatureTupleList;
}

// Function to assist in training and testing parent mass correction.  Given a spectrum,
// compute its self-convolution features, and return them as a list of tuples.
PyObject* PySpectrumGetPMCFeatures(PySpectrum* self, PyObject* args, PyObject* kwargs)
{
    PMCSpectrumInfo* SpectrumInfo;
    PMCInfo* Info;
    int PMCCount;
    PyObject* FeatureTupleList;
    PyObject* FeatureTuple;
    int FeatureIndex;
    int PMCFeatureCount = 64;
    //

    SpectrumInfo = GetPMCSpectrumInfo(self->Spectrum);
    ComputePMCFeatures(SpectrumInfo);

    // Count the PMCInfo nodes:
    PMCCount = 0;
    for (Info = SpectrumInfo->Head; Info; Info = Info->Next)
    {
        PMCCount++;
    }
    // Return a list of features.
    FeatureTupleList = PyList_New(0);
    for (Info = SpectrumInfo->Head; Info; Info = Info->Next)
    {
        FeatureTuple = PyTuple_New(PMCFeatureCount + 1);
        PyTuple_SetItem(FeatureTuple, 0, PyFloat_FromDouble(Info->ParentMass / (float)MASS_SCALE));
        for (FeatureIndex = 0; FeatureIndex < PMCFeatureCount; FeatureIndex++)
        {
            PyTuple_SetItem(FeatureTuple, FeatureIndex + 1, PyFloat_FromDouble(Info->Features[FeatureIndex]));
        }
        PyList_Append(FeatureTupleList, FeatureTuple);
        Py_DECREF(FeatureTuple); // Important!
    }
    FreePMCSpectrumInfo(SpectrumInfo);
    return FeatureTupleList;
}

// Function to assist in training and testing charge correction.  Given a spectrum,
// computes and return the charge-correction features.
PyObject* PySpectrumGetCCFeatures(PySpectrum* self, PyObject* args, PyObject* kwargs)
{
    int Charge1Flag = 0;
    PMCSpectrumInfo* SpectrumInfo1;
    PMCSpectrumInfo* SpectrumInfo2;
    PMCSpectrumInfo* SpectrumInfo3;
    float CCFeatures[64];
    int FeatureIndex;
    PyObject* FeatureTuple;
    //
    if (!PyArg_ParseTuple(args, "i", &Charge1Flag))
    {
        return NULL;
    }
    memset(CCFeatures, 0, sizeof(float) * 64);
    /////////////////////////////////
    // Charge 1 PMC:
    self->Spectrum->Charge = 1;
    self->Spectrum->ParentMass = (self->Spectrum->MZ * 1);
    SpectrumInfo1 = GetPMCSpectrumInfo(self->Spectrum);
    PerformPMC(SpectrumInfo1);
    /////////////////////////////////
    // Charge 2 PMC:
    self->Spectrum->Charge = 2;
    self->Spectrum->ParentMass = (self->Spectrum->MZ * 2) - HYDROGEN_MASS;
    SpectrumInfo2 = GetPMCSpectrumInfo(self->Spectrum);
    PerformPMC(SpectrumInfo2);
    /////////////////////////////////
    // Charge 3 PMC:
    self->Spectrum->Charge = 3;
    self->Spectrum->ParentMass = (self->Spectrum->MZ * 3) - 2 * HYDROGEN_MASS;
    SpectrumInfo3 = GetPMCSpectrumInfo(self->Spectrum);
    PerformPMC(SpectrumInfo3);
    if (Charge1Flag == 1)
    {
        //////////////////////////////////
        // Get features:
        GetChargeCorrectionFeatures1(SpectrumInfo1, SpectrumInfo2, SpectrumInfo3, CCFeatures);
    }
    else
    {
        GetChargeCorrectionFeatures2(SpectrumInfo2, SpectrumInfo3, CCFeatures);
    }
    FeatureTuple = PyTuple_New(64);
    for (FeatureIndex = 0; FeatureIndex < 64; FeatureIndex++)
    {
        PyTuple_SetItem(FeatureTuple, FeatureIndex, PyFloat_FromDouble(CCFeatures[FeatureIndex]));
    }

    FreePMCSpectrumInfo(SpectrumInfo1);
    FreePMCSpectrumInfo(SpectrumInfo2);
    FreePMCSpectrumInfo(SpectrumInfo3);
    return FeatureTuple;
}


PyObject* PySpectrumSetCharge(PySpectrum* self, PyObject* args, PyObject* kwargs)
{
    int NewCharge;
    if (!PyArg_ParseTuple(args, "i", &NewCharge))
    {
        return NULL;
    }
    self->Spectrum->Charge = NewCharge;
    // Reset the parent mass, based upon the M/Z:
    self->Spectrum->ParentMass = (self->Spectrum->MZ * NewCharge) - (NewCharge - 1) * HYDROGEN_MASS;
    // Set tweaks:
    TweakSpectrum(self->Spectrum->Node); //sam, comment out for phos pmc
    Py_INCREF(Py_None);
    return Py_None;
}

PyObject* PySpectrumGetMZ(PySpectrum* self, PyObject* args, PyObject* kwargs)
{
    float MZ;
    //
    MZ = self->Spectrum->MZ / (float)MASS_SCALE;
    return PyFloat_FromDouble(MZ);
}

PyObject* PySpectrumGetParentMass(PySpectrum* self, PyObject* args, PyObject* kwargs)
{
    float Mass;
    //
    Mass = self->Spectrum->ParentMass / (float)MASS_SCALE;
    return PyFloat_FromDouble(Mass);
}

PyObject* PySpectrumSetParentMass(PySpectrum* self, PyObject* args, PyObject* kwargs)
{
    float NewParentMass;
    int Charge = 0;
    MSSpectrum* Spectrum = self->Spectrum;
    //
    if (!PyArg_ParseTuple(args, "f|i", &NewParentMass, &Charge))
    {
        return NULL;
    }
    Spectrum->ParentMass = (int)(NewParentMass * MASS_SCALE + 0.5);
    if (Charge)
    {
        Spectrum->Charge = Charge;
    }
    if (Spectrum->Charge)
    {
        Spectrum->MZ = (Spectrum->ParentMass + (Spectrum->Charge - 1) * HYDROGEN_MASS) / Spectrum->Charge;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

PyObject* PySpectrumBYConvolve(PySpectrum* self, PyObject* args, PyObject* kwargs)
{
    // If TriplyChargedFlag is true, then we're looking for pairs with one singly-charged
    // peak and one doubly-charged peak.  Otherwise we're looking for two singly-charged peaks.
    int TriplyChargedFlag = 0;
    int Offset;
    float FloatOffset = 1.0;
    float Convolution = 0;
    int PeakIndex;
    int OtherMass;
    int Bin;
    float Intensity;
    //
    MSSpectrum* Spectrum = self->Spectrum;
    //
    if (!PyArg_ParseTuple(args, "|fi", &FloatOffset, &TriplyChargedFlag))
    {
        return NULL;
    }

    // Special case: If TriplyChargedFlag is -1, then compute direct self-convolution!
    // The direct self-convolution (dot product with self) is useful for scaling the
    // b,y convolutions.
    if (TriplyChargedFlag < 0)
    {
        Convolution = 0;
        for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
        {
            Bin = (Spectrum->Peaks[PeakIndex].Mass + 50) / 100;
            if (Bin >= 0 && Bin < Spectrum->IntensityBinCount)
            {
                Intensity = Spectrum->BinnedIntensitiesTight[Bin];
                Convolution += Spectrum->Peaks[PeakIndex].Intensity * Intensity; // * PeakScalingFactor;
            }
        }
        return PyFloat_FromDouble(Convolution);
    }
    Offset = (int)(FloatOffset * MASS_SCALE + 0.5);

    for (PeakIndex = 0; PeakIndex < Spectrum->PeakCount; PeakIndex++)
    {
        if (TriplyChargedFlag)
        {
            OtherMass = Spectrum->ParentMass + 2 * HYDROGEN_MASS - (2 * Spectrum->Peaks[PeakIndex].Mass) + Offset;
        }
        else
        {
            OtherMass = Spectrum->ParentMass + HYDROGEN_MASS - Spectrum->Peaks[PeakIndex].Mass + Offset;
        }
        Bin = ((OtherMass + 50) / 100);
        if (Bin < 0 || Bin >= Spectrum->IntensityBinCount)
        {
            continue;
        }
        Convolution += Spectrum->Peaks[PeakIndex].Intensity * Spectrum->BinnedIntensitiesTight[Bin];
    }
    return PyFloat_FromDouble(Convolution);
}

PyObject* PySpectrumCorrectCharge(PySpectrum* self, PyObject* args, PyObject* kwargs)
{
    // If TriplyChargedFlag is true, then we're looking for pairs with one singly-charged
    // peak and one doubly-charged peak.  Otherwise we're looking for two singly-charged peaks.
    int TriplyChargedFlag = 0;
    float FloatOffset = 1.0;
    float Convolution = 0;
    int Result;
    int ReturnScoresFlag = 0;
    float Model1Score;
    float Model2Score;
    //
    MSSpectrum* Spectrum = self->Spectrum;
    //
    if (!PyArg_ParseTuple(args, "|i", &ReturnScoresFlag))
    {
        return NULL;
    }
    Result = ChargeCorrectSpectrum(self->Spectrum->Node, &Model1Score, &Model2Score);
    if (ReturnScoresFlag)
    {
        return Py_BuildValue("ff", Model1Score, Model2Score);
    }
    else
    {
        return PyInt_FromLong(Result);
    }
}

void PySpectrumReportTagsFromTrie(PyObject* TagList, TrieNode* Root)
{
    TrieTagHanger* Hanger;
    TrieTag* Tag;
    int ChildIndex;
    TrieNode* Node;
    PyObject* TagTuple;
    //
    for (Hanger = Root->FirstTag; Hanger; Hanger = Hanger->Next)
    {
        Tag = Hanger->Tag;
        TagTuple = Py_BuildValue("fsf", Tag->PrefixMass / (float)DALTON,
            Tag->Tag, Tag->SuffixMass / (float)DALTON);
        PyList_Append(TagList, TagTuple);
    }

    for (ChildIndex = 0; ChildIndex < AMINO_ACIDS; ChildIndex++)
    {
        if (ChildIndex == 'I' - 'A' || ChildIndex == 'Q' - 'A')
        {
            continue; 
        }
        Node = Root->Children[ChildIndex];
        if (Node)
        {
            PySpectrumReportTagsFromTrie(TagList, Node);
        }
    }
}

int WriteTagsToList(TrieNode* Root, TrieTag* Tags, int MaxCount, int CurrentCount)
{
    //int RunningTotal;
    int ChildIndex;
    TrieNode* Child;
    TrieTagHanger* Hanger;
    //
    if (CurrentCount >= MaxCount)
    {
        return CurrentCount;
    }

    for (ChildIndex = 0; ChildIndex < AMINO_ACIDS; ChildIndex++)
    {
        if (ChildIndex == 'I' - 'A' || ChildIndex == 'Q' - 'A')
        {
            continue; 
        }

        Child = Root->Children[ChildIndex];
        if (Child)
        {
            CurrentCount = WriteTagsToList(Child, Tags, MaxCount, CurrentCount);
        }
    }
    for (Hanger = Root->FirstTag; Hanger; Hanger = Hanger->Next)
    {
        memcpy(Tags + CurrentCount, Hanger->Tag, sizeof(TrieTag));
        CurrentCount++;
        if (CurrentCount >= MaxCount)
        {
            return CurrentCount;
        }
    }
    return CurrentCount;
}

// Generate tags for this spectrum.  Optionally, use the new PRM scoring model
// to do so.
PyObject* PySpectrumGenerateTags(PySpectrum* self, PyObject* args, PyObject* kwargs)
{
    int CustomScoringModelFlag = 0;
    int Charge;
    MSSpectrum* Spectrum;
    SpectrumTweak* Tweak;
    TrieTag* Tags;
    PyObject* TagList;
    PyObject* TagTuple;
    int TagCount;
    PRMBayesianModel* Model;
    int MaximumTagCount = 200;
    int TagIndex;
    TrieTag* Tag;
    TrieNode* Root = NULL;
    float TagEdgeScoreMultiplier = 1.0;
    //static TrieTag* SortedFilteredTagList = NULL;
    //
    if (!PyArg_ParseTuple(args, "i|if", &Charge, &CustomScoringModelFlag, &TagEdgeScoreMultiplier))
    {
        return NULL;
    }
    Spectrum = self->Spectrum;
    Spectrum->Charge = Charge;
    TweakSpectrum(Spectrum->Node);
    Tweak = Spectrum->Node->Tweaks + (Charge - 1) * 2;
    if (Spectrum->Graph)
    {
        FreeTagGraph(Spectrum->Graph);
    }
    Spectrum->Graph = ConstructTagGraph(Spectrum);
    TagGraphAddNodes(Spectrum->Graph, Spectrum);
    // Look up the correct PRM scoring model:
    if (CustomScoringModelFlag)
    {
        Model = InteractiveBN;
    }
    else
    {
        // Use the current production model to score the nodes:
        if (Spectrum->Charge > 2)
        {
            Model = TAGModelCharge3;
        }
        else
        {
            Model = TAGModelCharge2;
        }
    }
    PrepareSpectrumForIonScoring(Model, Spectrum, 1);
    // Score PRMs using this model:
    TagGraphScorePRMNodes(Model, Spectrum->Graph, Spectrum, Tweak);
    
    TagGraphPopulateEdges(Spectrum->Graph);
    Tags = TagGraphGenerateTags(Spectrum->Graph, Spectrum, &TagCount, MaximumTagCount, Tweak, TagEdgeScoreMultiplier, InteractiveBN);
    
    // Note: This array of tags may have some duplicates.  Let's just build a trie and
    // use THAT for our output!
    Root = BuildTrieFromTags(Tags, TagCount, Root, MaximumTagCount);
    //if (!SortedFilteredTagList)
    //{
    //    SortedFilteredTagList = (TrieTag*)calloc(500, sizeof(TrieTag));
    //}
    TagCount = WriteTagsToList(Root, Tags, 500, 0);
    qsort(Tags, TagCount, sizeof(TrieTag), (QSortCompare)CompareTagScores); 
    TagList = PyList_New(0);
    for (TagIndex = 0; TagIndex < TagCount; TagIndex++)
    {
        Tag = Tags + TagIndex;
        TagTuple = Py_BuildValue("fsffii", Tag->PrefixMass / (float)DALTON,
            Tag->Tag, Tag->SuffixMass / (float)DALTON, Tag->Score, 
            Tag->TotalSkew, Tag->TotalAbsSkew);
        PyList_Append(TagList, TagTuple);
    }

    //PySpectrumReportTagsFromTrie(TagList, Root);
    FreeTrieNode(Root);
    
    return TagList;
}

PyObject* PySpectrumGetPRMScore(PySpectrum* self, PyObject* args, PyObject* kwargs)
{
    float Mass;
    int CustomScoringModelFlag = 0;
    PRMBayesianModel* Model;
    SpectrumTweak* Tweak;
    SpectrumNode* SpecNode;
    int IntMass;
    float Score;
    int VerboseFlag = 0;
    //
    if (!PyArg_ParseTuple(args, "f|ii", &Mass, &CustomScoringModelFlag, &VerboseFlag))
    {
        return NULL;
    }
    if (CustomScoringModelFlag)
    {
        Model = InteractiveBN;
    }
    else
    {
        if (self->Spectrum->Charge < 3)
        {
            Model = PRMModelCharge2;
        }
        else
        {
            Model = PRMModelCharge3;
        }
    }
    PrepareSpectrumForIonScoring(Model, self->Spectrum, 1);
    IntMass = (int)(Mass * MASS_SCALE);
    SpecNode = self->Spectrum->Node;
    if (SpecNode->Tweaks[2].Charge)
    {
        Tweak = SpecNode->Tweaks + 2;
    }
    else if (SpecNode->Tweaks[4].Charge)
    {
        Tweak = SpecNode->Tweaks + 4;
    }
    else
    {
        Tweak = SpecNode->Tweaks;
    }
    Score = GetIonPRMFeatures(self->Spectrum, Tweak, Model, IntMass, VerboseFlag);
    return PyFloat_FromDouble(Score);
    //GetIonPRMFeatures(self->Spectrum, Tweak, Model, Mass, 1);
}

PyObject* PySpectrumPlotPRMScores(PySpectrum* self, PyObject* args, PyObject* kwargs)
{
    char* FileName;
    FILE* PlotFile;
    int PRM;
    float Score;
    SpectrumTweak* Tweak;
    PRMBayesianModel* Model;
    int UseCustomModelFlag = 0;
    //
    if (!PyArg_ParseTuple(args, "s|i", &FileName, &UseCustomModelFlag))
    {
        return NULL;
    }
    if (UseCustomModelFlag)
    {
        Model = InteractiveBN; 
    }
    else
    {
        if (self->Spectrum->Charge < 3)
        {
            Model = PRMModelCharge2;
        }
        else
        {
            Model = PRMModelCharge3;
        }
    }

    PlotFile = fopen(FileName, "wb");
    if (!PlotFile)
    {
        sprintf(PythonErrorString, "Unable to open '%s'", FileName);
        ReportPythonError();
        return NULL;
    }
    
    PrepareSpectrumForIonScoring(Model, self->Spectrum, 1);
    Tweak = self->Spectrum->Node->Tweaks + (self->Spectrum->Charge * 2) - 2;
    for (PRM = 0; PRM < self->Spectrum->ParentMass; PRM += (DALTON / 10))
    {
        Score = GetIonPRMFeatures(self->Spectrum, Tweak, Model, PRM, 0);
        fprintf(PlotFile, "%.2f\t%.2f\t\n", PRM / (float)DALTON, Score);
    }
    fclose(PlotFile);
    Py_INCREF(Py_None);
    return Py_None;
}

void VerboseReportTopTag(TrieTag* Tag, Peptide* Match, MSSpectrum* Spectrum)
{
    int PRM;
    int AminoIndex;
    float PRMScore;
    int ValidFlag;
    int TagNodeIndex;
    TagGraphNode* Node;
    int Diff;
    int PeptidePRM;
    int AminoCount;
    int ModIndex;

    PRM = Tag->PrefixMass;
    AminoCount = strlen(Match->Bases);
    printf("  Tag %s %.2f (prefix %.2f, suffix %.2f)\n", Tag->Tag, Tag->Score, Tag->PrefixMass / (float)DALTON, Tag->SuffixMass / (float)DALTON);
    for (TagNodeIndex = 0; TagNodeIndex < 4; TagNodeIndex++)
    {
        // PRM is our node's mass.  
        // First question: Is it correct?
        ValidFlag = 0;
        PeptidePRM = 0;
        for (AminoIndex = 0; AminoIndex < AminoCount; AminoIndex++)
        {
            Diff = abs(PeptidePRM - PRM);
            if (Diff < GlobalOptions->Epsilon)
            {
                ValidFlag = 1;
                break;
            }
            
            PeptidePRM += PeptideMass[Match->Bases[AminoIndex]];
            for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
            {
                if (Match->AminoIndex[ModIndex] == AminoIndex)
                {
                    PeptidePRM += Match->ModType[ModIndex]->RealDelta;
                }
            }
        }
        // Second question: What's the score of the node?
        Node = Tag->Nodes[TagNodeIndex];
        PRMScore = Node->Score;
        if (ValidFlag)
        {
            printf("    [RIGHT]");
        }
        else
        {
            printf("    [wrong]");
        }
        printf(" %.2f %.2f", PRM / (float)DALTON, PRMScore);
        if (TagNodeIndex < 3)
        {
            printf(" -%c- ", Tag->Tag[TagNodeIndex]);
        }
        printf("\n");
        // Increment the PRM.
        PRM += PeptideMass[Tag->Tag[TagNodeIndex]];
        
    }
}

void VerboseReportTrueTagPRMs(Peptide* Match, MSSpectrum* Spectrum)
{
    int TruePRM;
    int PrevPRM;
    int AminoIndex;
    int AminoCount;
    int ModIndex;
    int Diff;
    TagGraphNode* BackNode;
    TagGraphNode* Node;
    TagGraphNode* BestNode;
    TagGraphNode* FirstNode;
    TagGraphNode* OldFirstNode = NULL;
    TagGraphEdge* BestEdge;
    float BestEdgeScore;
    TagGraphEdge* Edge;
    int BestDiff;
    float Score;
    //
    TruePRM = 0;
    AminoCount = strlen(Match->Bases);
    for (AminoIndex = 0; AminoIndex <= AminoCount; AminoIndex++)
    {
        BestNode = NULL;
        FirstNode = NULL;

        // Look for a close node:
        for (Node = Spectrum->Graph->FirstNode; Node; Node = Node->Next)
        {
            Diff = abs(Node->Mass - TruePRM);
            if (Diff < GlobalOptions->Epsilon)
            {
                if (!FirstNode)
                {
                    FirstNode = Node;
                }
                if (!BestNode || Node->Score > BestNode->Score)
                {
                    BestNode = Node;
                    BestDiff = Diff;
                }
            }
        }

        // Also, is there an edge to the PREVIOUS node?
        if (TruePRM)
        {
            BestEdge = NULL;
            BestEdgeScore = -9999;
            for (BackNode = OldFirstNode; BackNode; BackNode = BackNode->Next)
            {
                if (BackNode->Mass > PrevPRM + GlobalOptions->Epsilon)
                {
                    break;
                }
                for (Edge = BackNode->FirstEdge; Edge; Edge = Edge->Next)
                {
                    Diff = abs(Edge->ToNode->Mass - TruePRM);
                    if (Diff < GlobalOptions->Epsilon)
                    {
                        Score = Edge->FromNode->Score + Edge->ToNode->Score;
                        if (Score > BestEdgeScore)
                        {
                            BestEdge = Edge;
                            BestEdgeScore = Score;
                        }
                    }
                }
            }
            if (BestEdge)
            {
                printf("-Edge: From %.2f (%.2f) to %.2f (%.2f) via %c (skew %.2f)\n", BestEdge->FromNode->Mass / (float)DALTON,
                    BestEdge->FromNode->Score, BestEdge->ToNode->Mass / (float)DALTON, BestEdge->ToNode->Score,
                    BestEdge->Jump->Amino, BestEdge->Skew / (float)DALTON);
            }
            else
            {
                printf("-(no edge)\n");
            }
        }
        if (BestNode)
        {
            printf("PRM %.2f: Best node score PRM %.2f (diff %.2f) score %.2f\n",
                TruePRM / (float)DALTON, BestNode->Mass / (float)DALTON, BestDiff / (float)DALTON, BestNode->Score);
        }
        else
        {
            printf("PRM %.2f: (no node)\n", TruePRM / (float)DALTON);
        }
        // Add mass for this aa:
        OldFirstNode = FirstNode;
        PrevPRM = TruePRM;
        TruePRM += PeptideMass[Match->Bases[AminoIndex]];
        for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
        {
            if (Match->AminoIndex[ModIndex] == AminoIndex)
            {
                TruePRM += Match->ModType[ModIndex]->RealDelta;
            }
        }
        
    }
}

PyObject* PySpectrumCheckTagging(PySpectrum* self, PyObject* args, PyObject* kwargs)
{
    char* Annotation;
    Peptide* Match;
    MSSpectrum* Spectrum;
    PRMBayesianModel* Model;
    SpectrumTweak* Tweak;
    int Charge;
    int TagCount;
    int MaximumTagCount = 500;
    TrieTag* Tags;
    float TagEdgeScoreMultiplier = 10.0;
    int TagIndex;
    //
    if (!PyArg_ParseTuple(args, "s", &Annotation))
    {
        return NULL;
    }
    Match = GetPeptideFromAnnotation(Annotation);
    // Bail out if the match isn't valid:
    if (!Match)
    {
        Py_INCREF(Py_None);
        return Py_None;
    }
    Spectrum = self->Spectrum;
    Charge = self->Spectrum->Charge;
    // Look up the correct PRM scoring model:
    Model = InteractiveBN;
    PrepareSpectrumForIonScoring(Model, Spectrum, 1);
    Tweak = Spectrum->Node->Tweaks + (Charge - 1) * 2;
    // Generate some tags:
    if (Spectrum->Graph)
    {
        FreeTagGraph(Spectrum->Graph);
    }
    Spectrum->Graph = ConstructTagGraph(Spectrum);
    TagGraphAddNodes(Spectrum->Graph, Spectrum);
    // Score PRMs using this model:
    TagGraphScorePRMNodes(Model, Spectrum->Graph, Spectrum, Tweak);
    TagGraphPopulateEdges(Spectrum->Graph);
    Tags = TagGraphGenerateTags(Spectrum->Graph, Spectrum, &TagCount, MaximumTagCount, Tweak, TagEdgeScoreMultiplier, Model);
    printf("\nCheck tagging: Peptide %s\n", Annotation);
    printf("%s:%d\n", self->FileName, Spectrum->Node->FilePosition);
    ///////////////////////////////////////////////////////////////////////
    // Report on the tag graph node scores which correpsond to the true PRMs.  Are 
    // we missing a PRM entirely?  Missing an edge entirely?  Did we miss a tag
    // simply because the scores were medicore?
    VerboseReportTrueTagPRMs(Match, Spectrum);

    ///////////////////////////////////////////////////////////////////////
    // Report the top 10 tags.  Are they correct?  Partially correct?
    //printf("Top 10 tags (true peptide %s):\n", Annotation);
    for (TagIndex = 0; TagIndex < min(10, TagCount); TagIndex++)
    {
        VerboseReportTopTag(Tags + TagIndex, Match, Spectrum);
    }
    
    // Cleanup:
    FreePeptideNode(Match);
    Py_INCREF(Py_None);
    return Py_None;
}

PyObject* PySpectrumGetCutScores(PySpectrum* self, PyObject* args, PyObject* kwargs)
{
    Peptide* Match;
    char* Annotation;
    MSSpectrum* Spectrum;
    int NodeIndex;
    int AminoIndex;
    int PRM;
    float Score;
    int AminoCount;
    int ModIndex;
    PyObject* ReturnList;
    PRMBayesianNode* Node;
    //
    if (!PyArg_ParseTuple(args, "s", &Annotation))
    {
        return NULL;
    }
    Spectrum = self->Spectrum;
    Match = GetPeptideFromAnnotation(Annotation);
    // Force the spectrum's parent mass to match the right parent mass:
    Spectrum->ParentMass = Match->ParentMass;
    PrepareSpectrumForIonScoring(InteractiveBN, Spectrum, 1);
    AminoCount = strlen(Match->Bases);
    for (NodeIndex = 0, Node = InteractiveBN->Head; Node; NodeIndex++, Node = Node->Next)
    {
        PRM = 0;
        for (AminoIndex = 0; AminoIndex <= AminoCount; AminoIndex++)
        {
            ///////////////////////////////////////////////////////////////////////////////////////
            // Set values, and accumulate table entries:
            Node->Values[AminoIndex] = IonScoringGetNodeValue(InteractiveBN, Node, Spectrum, PRM, Match, AminoIndex);
            ///////////////////////////////////////////////////////////////////////////////////////
            // Add to PRM:
            if (AminoIndex == AminoCount)
            {
                break;
            }
            PRM += PeptideMass[Match->Bases[AminoIndex]];
            for (ModIndex = 0; ModIndex < MAX_PT_MODS; ModIndex++)
            {
                if (Match->AminoIndex[ModIndex] == AminoIndex)
                {
                    PRM += Match->ModType[ModIndex]->RealDelta;
                }
            }
        } // Amino loop
    } // NodeIndex loop

    // Iterate over the values arrays, building the return-list.
    ReturnList = PyList_New(0);
    for (AminoIndex = 0; AminoIndex <= AminoCount; AminoIndex++)
    {
        Score = PRMBNGetCutScore(Spectrum, InteractiveBN, AminoIndex);
        PyList_Append(ReturnList, PyFloat_FromDouble(Score));
    }

    // Cleanup:
    FreePeptideNode(Match);

    // Return:
    return ReturnList;
}

PyObject* PySpectrumGetMatchFeatures(PySpectrum* self, PyObject* args, PyObject* kwargs)
{
    PyObject* ReturnList;
    char* Annotation;
    Peptide* Match;
    MSSpectrum* Spectrum;
    float MatchFeatures[256];
    int FeatureIndex;
    int FeatureCount;
    //
    if (!PyArg_ParseTuple(args, "s", &Annotation))
    {
        return NULL;
    }
    Spectrum = self->Spectrum;
    Match = GetPeptideFromAnnotation(Annotation);
    FeatureCount = GetPeptideMatchFeaturesFull(Spectrum, Match, MatchFeatures);
    FreePeptideNode(Match);
    ReturnList = PyList_New(0);
    for (FeatureIndex = 0; FeatureIndex < FeatureCount; FeatureIndex++)
    {
        PyList_Append(ReturnList, PyFloat_FromDouble(MatchFeatures[FeatureIndex]));
    }
    
    return ReturnList;
}

PyObject* PySpectrumPrepareIonScoring(PySpectrum* self, PyObject* args, PyObject* kwargs)
{
    int CustomModelFlag = 0;
    PRMBayesianModel* Model;
    MSSpectrum* Spectrum = self->Spectrum;
    if (!PyArg_ParseTuple(args, "i", &CustomModelFlag))
    {
        return NULL;
    }
    if (CustomModelFlag)
    {
        Model = InteractiveBN;
    }
    else
    {
        if (self->Spectrum->Charge < 3)
        {
            Model = PRMModelCharge2;
        }
        else
        {
            Model = PRMModelCharge3;
        }
    }
    PrepareSpectrumForIonScoring(Model, Spectrum, 1);
    Py_INCREF(Py_None);
    return Py_None;
}

