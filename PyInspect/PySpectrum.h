#ifndef PY_SPECTRUM_H
#define PY_SPECTRUM_H
// PySpectrum: Python wrapper for a Spectrum object.
#include "Python.h"
#include "structmember.h"
#include "Utils.h"
#include "Spectrum.h"

typedef struct
{
    PyObject_HEAD
    MSSpectrum* Spectrum;
    char FileName[MAX_FILENAME_LEN];
    PyObject* MatchList; // list of PyPeptide instances
    int PrevMass;
    //struct Peptide* FirstMatch; // list of Peptide instances for matches
    //struct Peptide* LastMatch;
} PySpectrum;

extern PyTypeObject PySpectrumType;

void SpectrumSetCharge(MSSpectrum* Spectrum, int Charge);
PyObject* PySpectrumGetPeakCount(PySpectrum* self, PyObject* args, PyObject* kwargs);
PyObject* PySpectrumInit(PySpectrum* self, PyObject *args, PyObject *kwds);
PyObject* PySpectrumNew(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject* PySpectrumScore(PySpectrum* self, PyObject* args, PyObject* kwargs);
PyObject* PySpectrumScoreDetailed(PySpectrum* self, PyObject* args, PyObject* kwargs);
PyObject* PySpectrumLabelPeaks(PySpectrum* self, PyObject* args, PyObject* kwargs);
PyObject* PySpectrumCorrectParentMass(PySpectrum* self, PyObject* args, PyObject* kwargs);
PyObject* PySpectrumGetPMCFeatures(PySpectrum* self, PyObject* args, PyObject* kwargs);
PyObject* PySpectrumSetCharge(PySpectrum* self, PyObject* args, PyObject* kwargs);
PyObject* PySpectrumBYConvolve(PySpectrum* self, PyObject* args, PyObject* kwargs);
PyObject* PySpectrumSetParentMass(PySpectrum* self, PyObject* args, PyObject* kwargs);
PyObject* PySpectrumGetParentMass(PySpectrum* self, PyObject* args, PyObject* kwargs);
PyObject* PySpectrumGetCCFeatures(PySpectrum* self, PyObject* args, PyObject* kwargs);
PyObject* PySpectrumCorrectCharge(PySpectrum* self, PyObject* args, PyObject* kwargs);
PyObject* PySpectrumGenerateTags(PySpectrum* self, PyObject* args, PyObject* kwargs);
PyObject* PySpectrumGetPRMScore(PySpectrum* self, PyObject* args, PyObject* kwargs);
PyObject* PySpectrumPlotPRMScores(PySpectrum* self, PyObject* args, PyObject* kwargs);
PyObject* PySpectrumCheckTagging(PySpectrum* self, PyObject* args, PyObject* kwargs);
PyObject* PySpectrumGetCutScores(PySpectrum* self, PyObject* args, PyObject* kwargs);
PyObject* PySpectrumGetMatchFeatures(PySpectrum* self, PyObject* args, PyObject* kwargs);
PyObject* PySpectrumPrepareIonScoring(PySpectrum* self, PyObject* args, PyObject* kwargs);
PyObject* PySpectrumGetMZ(PySpectrum* self, PyObject* args, PyObject* kwargs);

void PySpectrumDealloc(PyObject* selfobject);

// Methods of the PySpectrum class
static PyMethodDef PySpectrumMethods[] = 
{
    {"GetPeakCount", (PyCFunction)PySpectrumGetPeakCount, METH_VARARGS | METH_KEYWORDS, 
      "Return peak count for a scan"},
    {"ScorePeptide", (PyCFunction)PySpectrumScore, METH_VARARGS | METH_KEYWORDS, 
      "Score a peptide match for this spectrum"},
    {"ScorePeptideDetailed", (PyCFunction)PySpectrumScoreDetailed, METH_VARARGS | METH_KEYWORDS, 
      "Score a peptide match for this spectrum; return all the scoring features"},
    {"LabelPeaks", (PyCFunction)PySpectrumLabelPeaks, METH_VARARGS | METH_KEYWORDS, 
      "Label spectrum peaks using a peptide annotation"},
    {"CorrectParentMass", (PyCFunction)PySpectrumCorrectParentMass, METH_VARARGS | METH_KEYWORDS, 
      "Select correct charge and parent mass for the spectrum"},
    {"GetPMCFeatures", (PyCFunction)PySpectrumGetPMCFeatures, METH_VARARGS | METH_KEYWORDS, 
      "Compute parent-mass-correction features for the spectrum"},
    {"SetCharge", (PyCFunction)PySpectrumSetCharge, METH_VARARGS | METH_KEYWORDS, 
      "Adjust the spectrum's charge"},
    {"BYConvolve", (PyCFunction)PySpectrumBYConvolve, METH_VARARGS | METH_KEYWORDS, 
      "Perform b/y peak convolution"},
    {"SetParentMass", (PyCFunction)PySpectrumSetParentMass, METH_VARARGS | METH_KEYWORDS, 
      "Set the parent mass"},
    {"GetParentMass", (PyCFunction)PySpectrumGetParentMass, METH_VARARGS | METH_KEYWORDS, 
      "Returns the parent mass"},
    {"GetCCFeatures", (PyCFunction)PySpectrumGetCCFeatures, METH_VARARGS | METH_KEYWORDS, 
      "Compute charge correction features for the spectrum"},
    {"CorrectCharge", (PyCFunction)PySpectrumCorrectCharge, METH_VARARGS | METH_KEYWORDS, 
      "Get the corrected charge for this spectrum"},
    {"GenerateTags", (PyCFunction)PySpectrumGenerateTags, METH_VARARGS | METH_KEYWORDS, 
      "Generate tags for this spectrum"},
    {"GetPRMScore", (PyCFunction)PySpectrumGetPRMScore, METH_VARARGS | METH_KEYWORDS, 
      "Get the score for a prefix residue mass (PRM)"},
    {"PlotPRMScores", (PyCFunction)PySpectrumPlotPRMScores, METH_VARARGS | METH_KEYWORDS, 
      "Output a plot of PRM scores for this spectrum"},
    {"CheckTagging", (PyCFunction)PySpectrumCheckTagging, METH_VARARGS | METH_KEYWORDS, 
      "Test whether this spectrum can generate a tag for a specified peptide"},
    {"GetCutScores", (PyCFunction)PySpectrumGetCutScores, METH_VARARGS | METH_KEYWORDS, 
      "Compute cut-point scores for the specified peptide annotation"},
    {"GetMatchFeatures", (PyCFunction)PySpectrumGetMatchFeatures, METH_VARARGS | METH_KEYWORDS, 
      "Get features for scoring a peptide match"},
    {"PrepareIonScoring", (PyCFunction)PySpectrumPrepareIonScoring, METH_VARARGS | METH_KEYWORDS, 
      "Force a call to PrepareSpectrumForIonScoring"},
    {"GetMZ", (PyCFunction)PySpectrumGetMZ, METH_VARARGS | METH_KEYWORDS,
      "Get the m/z for this spectrum"},
      

    {NULL},
};


// Methods (currently none) of the MSPeak class
static PyMemberDef PySpectrumMembers[] =
{
    {NULL},
};

// Getters and setters for the PySpectrum class.  (Should be used
// if Python code will be modifying the run dynamically)
static PyGetSetDef PySpectrumGetSet[] = 
{
    {NULL}  // Sentinel 
};



#endif // PY_SPECTRUM_H
