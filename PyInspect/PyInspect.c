#include "CMemLeak.h"
#include "Python.h"
#include "Trie.h"
#include "Score.h"
#include "PySpectrum.h"
#include "BN.h"
#include "PyUtils.h"
#include "Errors.h"
#include "FreeMod.h"
#include "IonScoring.h"
#include "SVM.h"
#include "LDA.h"
#include "ParseInput.h"
#include "ParseXML.h"

PyObject* InspectError;
PRMBayesianModel* InteractiveBN = NULL;

PyObject* PyResetIonScoring(PyObject* self, PyObject* args)
{
    int IntensityScheme;
    float IntensityRadius;
    int CutFlag = 0;
    int NoiseModelFlag = 0;
    //
    if (!PyArg_ParseTuple(args, "if|ii", &IntensityScheme, &IntensityRadius, &CutFlag, &NoiseModelFlag))
    {
        return NULL;
    }
    FreePRMBayesianModel(InteractiveBN);
    InteractiveBN = (PRMBayesianModel*)calloc(1, sizeof(PRMBayesianModel));
    InteractiveBN->NoiseModel = NoiseModelFlag;
    InteractiveBN->IntensityScheme = IntensityScheme;
    switch (InteractiveBN->IntensityScheme)
    {
    case 0:
    case 1:
    case 4:
        InteractiveBN->MinIntensityLevel = 3; 
        break;
    case 2:
    case 3:
        InteractiveBN->MinIntensityLevel = 2; 
        break;
    default:
        REPORT_ERROR(0);
        break;
    }
    
    InteractiveBN->IntensityRadius = (int)(IntensityRadius * DALTON);
    InteractiveBN->HalfIntensityRadius = InteractiveBN->IntensityRadius / 2;
    InteractiveBN->CutFlag = CutFlag;
    Py_INCREF(Py_None);
    return Py_None;
}

PyObject* PyAddIonScoringNode(PyObject* self, PyObject* args)
{
    int NodeType;
    int NodeInfoA;
    char* Name;
    int FragmentType = evFragmentTypeNone;
    float NodeMassOffset;
    if (!PyArg_ParseTuple(args, "sii|fi", &Name, &NodeType, &NodeInfoA, &NodeMassOffset, &FragmentType))
    {
        return NULL;
    }
    // Add a node:
    AddPRMBayesianNode(InteractiveBN, Name, NodeType, NodeInfoA, NodeMassOffset, FragmentType);
    // Return the node's index:
    return PyInt_FromLong(InteractiveBN->NodeCount - 1);
}

PyObject* PySetIonScoringNodeParents(PyObject* self, PyObject* args)
{
    int NodeIndex;
    PyObject* ParentIndexList;
    PRMBayesianNode* Node;
    PRMBayesianNode* Parent;
    int OverallBlockSize;
    int ParentIndex;
    int OtherParentIndex;
    //
    if (!PyArg_ParseTuple(args, "iO", &NodeIndex, &ParentIndexList))
    {
        return NULL;
    }
    // Validate input:
    if (NodeIndex < 0 || NodeIndex >= InteractiveBN->NodeCount)
    {
        sprintf(PythonErrorString, "Illegal node index %d in SetIonScoringNodeParents", NodeIndex);
        ReportPythonError();
        return NULL;
    }
    Node = InteractiveBN->Nodes[NodeIndex];
    
    // Free the OLD parents, if any:
    SafeFree(Node->Parents);
    SafeFree(Node->ParentBlocks);

    // Set the parents of this node:
    Node->ParentCount = PyList_Size(ParentIndexList);
    if (Node->ParentCount)
    {
        Node->Parents = (PRMBayesianNode**)calloc(sizeof(PRMBayesianNode*), Node->ParentCount);
        Node->ParentBlocks = (int*)calloc(sizeof(int), Node->ParentCount);
    }
    for (ParentIndex = 0; ParentIndex < Node->ParentCount; ParentIndex++)
    {
        Parent = InteractiveBN->Nodes[PyInt_AsLong(PyList_GetItem(ParentIndexList, ParentIndex))];
        Node->Parents[ParentIndex] = Parent;
    }
    // Set the parent block sizes.  Node->ParentBlocks[n] is the total number of combinations
    // for values of parents n+1 and beyond (or 1, if n is ParentCount - 1).  When indexing into
    // the probability tables, we use these blocks.
    for (ParentIndex = 0; ParentIndex < Node->ParentCount; ParentIndex++)
    {
        OverallBlockSize = Node->ValueCount;
        for (OtherParentIndex = ParentIndex + 1; OtherParentIndex < Node->ParentCount; OtherParentIndex++)
        {
            OverallBlockSize *= Node->Parents[OtherParentIndex]->ValueCount;
        }
        Node->ParentBlocks[ParentIndex] = OverallBlockSize;
        Parent = Node->Parents[ParentIndex];
        OverallBlockSize /= Parent->ValueCount;
    }
    // Allocate probability tables:
    Node->TableSize = Node->ValueCount;
    if (Node->ParentCount)
    {
        Node->TableSize = (Node->Parents[0]->ValueCount * Node->ParentBlocks[0]);
    }
    SafeFree(Node->CountTable);
    Node->CountTable = (int*)calloc(Node->TableSize, sizeof(int));
    SafeFree(Node->ProbTable);
    Node->ProbTable = (float*)calloc(Node->TableSize, sizeof(float));
    Py_INCREF(Py_None);
    return Py_None;
}

void TrainNoiseModelRandomMasses(PRMBayesianModel* Model, MSSpectrum* Spectrum)
{
    int Bin;
    int BinIndex;
    //
    for (BinIndex = 0; BinIndex < 20; BinIndex++)
    {
        Bin = rand() % Spectrum->IntensityBinCount;
        Model->RandomIntensityCounts[Spectrum->BinnedIntensityLevels[Bin]]++;
    }
}

PyObject* PyTrainBNOnSpectrum(PyObject* self, PyObject* args)
{
    PySpectrum* SpectrumObject;
    char* PeptideAnnotation;
    Peptide* Match;
    int PRM = 0;
    int AminoCount;
    int AminoIndex;
    int ModIndex;
    int NodeIndex;
    int TableIndex;
    int ParentIndex;
    MSSpectrum* Spectrum;
    PRMBayesianNode* Node;
    //
    if (!PyArg_ParseTuple(args, "Os", &SpectrumObject, &PeptideAnnotation))
    {
        return NULL;
    }
    Match = GetPeptideFromAnnotation(PeptideAnnotation);
    if (!Match)
    {
        REPORT_ERROR(0);
        return NULL;
    }
    Spectrum = SpectrumObject->Spectrum;
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
            // Accumulate PRM from the prefix so far:
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

    // Iterate over the values arrays, accumulate counts in the frequency tables:
    for (NodeIndex = 0; NodeIndex < InteractiveBN->NodeCount; NodeIndex++)
    {
        for (AminoIndex = 0; AminoIndex <= AminoCount; AminoIndex++)
        {
            Node = InteractiveBN->Nodes[NodeIndex];
            TableIndex = 0;
            for (ParentIndex = 0; ParentIndex < Node->ParentCount; ParentIndex++)
            {
                TableIndex += Node->Parents[ParentIndex]->Values[AminoIndex] * Node->ParentBlocks[ParentIndex];
            }
            TableIndex += Node->Values[AminoIndex];
            if (TableIndex >= Node->TableSize)
            {
                // Panic!
                REPORT_ERROR(0);
                TableIndex = 0;
            }
            Node->CountTable[TableIndex]++;
        }
    }

    // And, count how frequent the various intensity levels are for a random mass:
    TrainNoiseModelRandomMasses(InteractiveBN, Spectrum);

    // Cleanup:
    FreePeptideNode(Match);

    // Return:
    Py_INCREF(Py_None);
    return Py_None;
}

PyObject* PyDebugPrintPRMBayesianModel(PyObject* self, PyObject* args)
{
    if (!PyArg_ParseTuple(args, ""))
    {
        return NULL;
    }

    DebugPrintPRMBayesianModel(InteractiveBN);
    Py_INCREF(Py_None);
    return Py_None;
}

// Save InteractiveBN to disk:
PyObject* PySaveBNModel(PyObject* self, PyObject* args)
{
    char* FileName;
    if (!PyArg_ParseTuple(args, "s", &FileName))
    {
        return NULL;
    }

    SavePRMBayesianModel(InteractiveBN, FileName);
    Py_INCREF(Py_None);
    return Py_None;
}

// Load InteractiveBN from file; return node count
PyObject* PyLoadBNModel(PyObject* self, PyObject* args)
{
    char* FileName;
    if (!PyArg_ParseTuple(args, "s", &FileName))
    {
        return NULL;
    }

    FreePRMBayesianModel(InteractiveBN);
    InteractiveBN = LoadPRMBayesianModel(FileName);
    if (!InteractiveBN)
    {
        return PyInt_FromLong(-1);
    }
    return PyInt_FromLong(InteractiveBN->NodeCount);
}

// Convert the COUNT tables of the bayesian network into PROBABILITY tables.
PyObject* PyComputeBNProbabilityTables(PyObject* self, PyObject* args)
{
    if (!PyArg_ParseTuple(args, ""))
    {
        return NULL;
    }
    ComputePRMBayesianModelProbabilityTables(InteractiveBN, 1);
    Py_INCREF(Py_None);
    return Py_None;
}

PyObject* PyComputeMutualInformation(PyObject* self, PyObject* args)
{
    PyObject* ReturnList;
    PyObject* NodeEntropyList;
    PRMBayesianNode* Node;
    PRMBayesianNode* Parent;
    float EntropySum[256];
    int TableIndex;
    float Entropy;
    float JointEntropy;
    int ParentIndex;
    int TempIndex;
    int ParentValue;
    float Probability;
    int Value;
    int ValueIndex;
    int NodeIndex;
    float NodeEntropy[512];
    float MutualInformation;
    int FullTableCount;
    if (!PyArg_ParseTuple(args, ""))
    {
        return NULL;
    }

    ReturnList = PyList_New(0);

    // Compute the entropy of each node:
    for (NodeIndex = 0, Node = InteractiveBN->Head; Node; Node = Node->Next, NodeIndex++)
    {
        ////////////////////////////////////////////////////////////////
        // Compute the node's entropy:
        memset(EntropySum, 0, sizeof(float) * 256);
        FullTableCount = 0;
        for (TableIndex = 0; TableIndex < Node->TableSize; TableIndex++)
        {
            FullTableCount += Node->CountTable[TableIndex];
        }
        for (TableIndex = 0; TableIndex < Node->TableSize; TableIndex++)
        {
            Probability = Node->CountTable[TableIndex] / (float)FullTableCount;
            Value = TableIndex % Node->ValueCount;
            EntropySum[Value] += Probability; 
        }
        Entropy = 0;
        for (ValueIndex = 0; ValueIndex < Node->ValueCount; ValueIndex++)
        {
            printf("Node %d %s value %d: Odds %.6f\n", NodeIndex, Node->Name, ValueIndex, EntropySum[ValueIndex]);
            if (EntropySum[ValueIndex] > 0.0)
            {
                Entropy -= EntropySum[ValueIndex] * (float)log(EntropySum[ValueIndex]);
            }
        }
        NodeEntropy[NodeIndex] = Entropy;
    }
    for (NodeIndex = 0, Node = InteractiveBN->Head; Node; Node = Node->Next, NodeIndex++)
    {
        NodeEntropyList = PyList_New(0);
        PyList_Append(ReturnList, NodeEntropyList);
        PyList_Append(NodeEntropyList, PyInt_FromLong(NodeIndex));
        PyList_Append(NodeEntropyList, PyString_FromString(Node->Name));
        PyList_Append(NodeEntropyList, PyFloat_FromDouble(NodeEntropy[NodeIndex]));
        ////////////////////////////////////////////////////////////////
        // Compute the node's joint entropy with each parent:
        for (ParentIndex = 0; ParentIndex < Node->ParentCount; ParentIndex++)
        {
            // EntropySum[ParentValue*Block+Value] = probability that parent has ParentValue
            // and node has Value
            Parent = Node->Parents[ParentIndex];
            memset(EntropySum, 0, sizeof(float) * 256);
            for (TableIndex = 0; TableIndex < Node->TableSize; TableIndex++)
            {
                TempIndex = TableIndex;
                if (ParentIndex)
                {
                    TempIndex = TempIndex % Node->ParentBlocks[ParentIndex - 1];
                }
                ParentValue = TempIndex / Node->ParentBlocks[ParentIndex];
                Probability = Node->CountTable[TableIndex] / (float)FullTableCount;
                //Probability = (float)exp(Node->ProbTable[TableIndex]);
                Value = TableIndex % Node->ValueCount;
                EntropySum[ParentValue * Node->ValueCount + Value] += Probability; 
            }
            JointEntropy = 0;
            for (ValueIndex = 0; ValueIndex < (Node->ValueCount * Parent->ValueCount); ValueIndex++)
            {
                ParentValue = ValueIndex / Node->ValueCount;
                printf("Node %d %s value %d and parent (%d %s) has value %d: Odds %.6f\n", NodeIndex, 
                    Node->Name, 
                    ValueIndex % Node->ValueCount, Parent->Index, Parent->Name, ParentValue, EntropySum[ValueIndex]);
                if (EntropySum[ValueIndex] > 0.0)
                {
                    JointEntropy -= EntropySum[ValueIndex] * (float)log(EntropySum[ValueIndex]);
                }
            }
            MutualInformation = (NodeEntropy[NodeIndex] + NodeEntropy[Parent->Index] - JointEntropy);
            printf("Node %d(%s) parent %d(%s):\n", NodeIndex, Node->Name, Parent->Index, Parent->Name);
            printf("  Child entropy %.6f, parent entropy %.6f\n", NodeEntropy[NodeIndex], NodeEntropy[Parent->Index]);
            printf("  Joint entropy: %.6f\n", JointEntropy);
            printf("  Mutual information: %.6f\n", MutualInformation);
            printf("  Conditional entropy (Child|Parent): %.6f\n", JointEntropy - NodeEntropy[Parent->Index]);
            printf("  Normalized MI: %.6f\n", MutualInformation / NodeEntropy[NodeIndex]);
            PyList_Append(NodeEntropyList, PyFloat_FromDouble(MutualInformation / NodeEntropy[NodeIndex]));
        }
        Py_DECREF(NodeEntropyList);
    }
    return ReturnList;
}

PyObject* PyGetBNFeatureNames(PyObject* self, PyObject* args)
{
    PyObject* ReturnList;
    PRMBayesianNode* Node;
    //
    if (!PyArg_ParseTuple(args, ""))
    {
        return NULL;
    }
    ReturnList = PyList_New(0);
    for (Node = InteractiveBN->Head; Node; Node = Node->Next)
    {
        PyList_Append(ReturnList, PyString_FromString(Node->Name));
    }
    return ReturnList;
}

PyObject* PyComputeBNValuesForSpectrum(PyObject* self, PyObject* args)
{
    PySpectrum* SpectrumObject;
    char* PeptideAnnotation;
    Peptide* Match;
    int PRM = 0;
    int AminoCount;
    int AminoIndex;
    int ModIndex;
    int NodeIndex;
    MSSpectrum* Spectrum;
    PRMBayesianNode* Node;
    PyObject* ReturnList;
    PyObject* NodeValueList;
    //
    if (!PyArg_ParseTuple(args, "Os", &SpectrumObject, &PeptideAnnotation))
    {
        return NULL;
    }
    Match = GetPeptideFromAnnotation(PeptideAnnotation);
    if (!Match)
    {
        REPORT_ERROR(0);
        return NULL;
    }
    Spectrum = SpectrumObject->Spectrum;
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
        NodeValueList = PyList_New(0);
        PyList_Append(ReturnList, NodeValueList);
        for (NodeIndex = 0; NodeIndex < InteractiveBN->NodeCount; NodeIndex++)
        {
            Node = InteractiveBN->Nodes[NodeIndex];
            PyList_Append(NodeValueList, PyInt_FromLong(Node->Values[AminoIndex]));
            //Py_DECREF(NodeValueList);
        }
    }

    // Cleanup:
    FreePeptideNode(Match);

    // Return:
    return ReturnList;
}

PyObject* PyFinishIonScoringNetwork(PyObject* self, PyObject* args)
{
    if (!PyArg_ParseTuple(args, ""))
    {
        return NULL;
    }

    // Perform any activities necessary to finalizing InteractiveBN:
    BuildModelFlankList(InteractiveBN);
    Py_INCREF(Py_None);
    return Py_None;
}

PyObject* PyReloadPMC(PyObject* self, PyObject* args)
{
    //
    if (!PyArg_ParseTuple(args, "|"))
    {
        return NULL;
    }

    // Reload parent mass correction and charge-correction models:

#ifdef PMC_USE_SVM
    LoadPMCSVM(1);
#else
    LoadPMCLDA(1);
#endif    

#ifdef CC_USE_SVM
    LoadCCModelSVM(1);
#else
    LoadCCModelLDA(1);
#endif
    Py_INCREF(Py_None);
    return Py_None;
}

static PyMethodDef PyInspectMethods[] = 
{
    {"ResetIonScoring", PyResetIonScoring, 1, "Reset the ion scoring network"},
    {"AddIonScoringNode", PyAddIonScoringNode, 1, "Add a node to the ion scoring network"},
    {"SetIonScoringNodeParents", PySetIonScoringNodeParents, 1, "Set the parent(s) of an ion scoring node"},
    {"FinishIonScoringNetwork", PyFinishIonScoringNetwork, 1, "Finalize an ion scoring network"},
    {"TrainBNOnSpectrum", PyTrainBNOnSpectrum, 1, "Accumulate counts for network nodes, given a spectrum and peptide"},
    {"DebugPrintBNModel", PyDebugPrintPRMBayesianModel, 1, "Debug print"},
    {"SaveBNModel", PySaveBNModel, 1, "Save model to a binary file"},
    {"LoadBNModel", PyLoadBNModel, 1, "Load from binary file (as written by SaveBNModel)"},
    {"ComputeBNProbabilityTables", PyComputeBNProbabilityTables, 1, "Compute probability tables for a BNModel"},
    {"ComputeBNValuesForSpectrum", PyComputeBNValuesForSpectrum, 1, "Compute values for nodes in the BNModel"},
    {"GetBNFeatureNames", PyGetBNFeatureNames, 1, "Return a list of names of nodes in the bayesian network"},
    {"ComputeMutualInformation", PyComputeMutualInformation, 1, "Compute MutualInformation for nodes and their parents"},
    {"ReloadPMC", PyReloadPMC, 1, "Reset PMC / CC models"},
    //{"erf", PyErrorFunction, METH_VARARGS, "return the error function erf(x)"},
    //{"GammaIncomplete", PyGammaIncomplete, METH_VARARGS, "return the incomplete gamma function g(a, x)"},
    //{"foo", ex_foo, 1, "foo() doc string"},
    {NULL, NULL}
};

// Cleanup, called by Python when unloading.  Deallocate memory:
void PyInspectCleanup(void)
{
    FreeMassDeltaByMass();
    FreeMassDeltas();
    FreeIsSubDecoration();
    //FreeTaggingModel();
    FreeJumpingHash();
    FreeSVMModels();
    FreeBayesianModels();
    FreeLDAModels();
    FreePRMBayesianModel(InteractiveBN);
    InteractiveBN = NULL;
    SafeFree(GlobalOptions);
    GlobalOptions = NULL;
    FreeCCModelSVM();
    FreeTagSkewScores();
    SafeFree(MassDeltaByIndex);
    MassDeltaByIndex = NULL;
    FreeMZXMLParseCursor();
    FreeMZDataParseCursor();

}

PyMODINIT_FUNC initPyInspect(void)
{
    PyObject* Module;
    ////////////////////
    PySpectrumType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&PySpectrumType) < 0)
    {
        return;
    }
    Module = Py_InitModule("PyInspect", PyInspectMethods);

    // Add the Error object:
    InspectError = PyErr_NewException("Inspect.error", NULL, NULL);
    Py_INCREF(InspectError);
    PyModule_AddObject(Module, "error", InspectError);
    InitErrors();

    // Add the Spectrum object:
    PyModule_AddObject(Module, "Spectrum", (PyObject *)&PySpectrumType);

    // Create an ion scoring network, for interactive use:
    InteractiveBN = (PRMBayesianModel*)calloc(1, sizeof(PRMBayesianModel));

    // Perform some standard loading here, like amino acid masses.
    AllocMassDeltaByIndex();
    InitOptions();
    sprintf(GlobalOptions->ResourceDir, ".%c", SEPARATOR);
    LoadPeptideMasses(NULL);
    PeptideMass['C'] += 57000; // ASSUMED: All cysteines carry the +57 modification.
    LoadMassDeltas(NULL, 0);
    InitBayesianModels();
    SetTagSkewScores();
    //LoadFlankingAminoEffects();
    //LoadCCModel();
#ifdef MQSCORE_USE_SVM
    InitPValueSVM();
#else
    InitPValueLDA();
#endif
    PopulateJumpingHash();
    
    // Set the blind-flag to TRUE so that modified peptides 
    // incur a score-penalty:
    GlobalOptions->RunMode |= RUN_MODE_BLIND;
    // Register our cleanup function to run at exit:
    Py_AtExit(PyInspectCleanup);
}
