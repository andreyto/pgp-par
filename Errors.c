//Copyright 2007, The Regents of the University of California
//All Rights Reserved
//
//Permission to use, copy, modify and distribute any part of this 
//program for educational, research and non-profit purposes, without fee, 
//and without a written agreement is hereby granted, provided that the 
//above copyright notice, this paragraph and the following three paragraphs 
//appear in all copies.
//
//Those desiring to incorporate this work into commercial 
//products or use for commercial purposes should contact the Technology 
//Transfer & Intellectual Property Services, University of California, 
//San Diego, 9500 Gilman Drive, Mail Code 0910, La Jolla, CA 92093-0910, 
//Ph: (858) 534-5815, FAX: (858) 534-7345, E-MAIL:invent@ucsd.edu.
//
//IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY 
//FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, 
//INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN 
//IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY 
//OF SUCH DAMAGE.
//
//THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE UNIVERSITY 
//OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, 
//ENHANCEMENTS, OR MODIFICATIONS.  THE UNIVERSITY OF CALIFORNIA MAKES NO 
//REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR 
//EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
//MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF 
//THE SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.
#include "CMemLeak.h"
#include <stdio.h>
#include "Inspect.h"
#include "Errors.h"

void AssertionFailed(char* Assertion, char* FileName, int LineNumber)
{
    printf("** ASSERTION FAILED line %d file '%s':\n   '%s'\n", LineNumber, FileName, Assertion);
}

typedef struct NumberedError
{
    int ID;
    char* Message;
} NumberedError;

NumberedError ErrorMessages[] = {
    {0, "Unhandled exception"},
    {1, "Out of memory"},
    {2, "Out of disk space"},
    {3, "Missing required file '%s'"},
    {4, "Internal assertion '%s'"},
    {5, "File '%s' not found"}, 
    {6, "Error in LoadBayesianModel: Bogus feature count %d"},
    {7, "Bogus-looking probability table size %d for feature %d"},
    {8, "Unable to open requested file '%s'"},
    {9, "Scan number range (%d...%d) includes no spectra!"},
    {10, "Only %d top-scoring matches for charge state; not recalibrating the p-value curve."},
    {11, "No spectra were specified to search."},
    {12, "No GFF files were specified as input to build an MS2DB file."},
    {13, "Ignoring unknown command '%s' from Inspect input file"},
    {14, "Syntax error on line %d of file %s"},
    {15, "No valid exons found in GFF files"},
    {16, "Linked exons %d...%d and %d...%d have incompatible reading frames"},
    {17, "Consecutive GFF exons %d...%d and %d...%d come from same gene, but can't be linked because they overlap"},
    {18, "Invalid command-line argument '%s'"},
    {19, "Command-line argument '%s' requires a parameter."},
    {20, "Invalid coordinates %d...%d on line %d of file %s"},
    {21, "Length-1 exon at %d is a codon center, but doesn't link in and out"},
    {22, "Unable to cover GFF gene '%s'"},
    {23, "Coverage of GFF gene '%s' is incomplete"},
    {24, "Unhandled instance: %d"},
    {25, "XML parse exception: '%s'"},
    {26, "Error linking exons %d and %d in gene '%s'"},
    {27, "XML line %d: Unexpected tag '%s'"},
    {28, "XML line %d: Unexpected attribute '%s' for XML tag '%s'"},
    {29, "Exon %d of gene '%s' is too long!"},
    {30, "Spectrum file '%s' has an abnormal extension; attempting to treat it as a .dta file"},
    {31, "Too many peaks in spectrum (scan %d, file %s); dropping extras!"},
    {32, "Illegal peak mass in scan %d of file %s"},
    {33, "Syntax error on line %d of input file %s"},
    {34, "Modifications specified, but no PTMs permitted in peptides.  Use 'mods,1' to permit modified peptides."},
    {35, "Too many PTMs specified in input file - ignoring '%s'"},
    {36, "Invalid modification type '%s'"},
    {37, "Illegal amino acid specificity '%s' for modification"},
    {38, "Invalid tag length '%d': Valid values are 1 through 6"},
    {39, "Input file parameter '%s' doesn't take a value"},
    {40, "Input file parameter '%s' requires a string value"},
    {41, "Input file parameter '%s' requires an integer value"},
    {42, "Invalid mass %dDa at position %d in spectrum file '%s'"},
    {43, "Invalid mass %dDa in spectrum file"},
    {44, "Invalid mass ppm %d - should be in the range 1...1000"},
    {45, "Peak for spectrum %s:%d is %dDa - possible corruption"},
    {46, "Invalid scoring model specified: Charge must be 2 or 3"},
    {47, "Invalid RequiredTermini value '%d' specified"},
    {48, "Peak for spectrum %s:%d has intensity %f - possible corruption"},
    {49, "Out of memory - failed to allocate %d bytes"},
    {-1, NULL}
};

int ErrorMessageCount;

void InitErrors()
{
    ErrorMessageCount = sizeof(ErrorMessages) / sizeof(char*);
}

// Report an error - write it to GlobalOptions->ErrorFile (if GlobalOptions exists)
// and to stderr, and increment the count of reported errors.
void ReportError(int ErrorSeverity, int ErrorID, int SourceLine, char* SourceFileName, int ArgType,
                 const char* StrA, const char* StrB, const char* StrC, 
                 int IntA, int IntB, int IntC, int IntD,
                 float FloatA, float FloatB)
{
    char* ErrorMessage;
    int ErrorIndex;
    FILE* ErrorFile;
    FILE* ErrorFile2;
    //
    if (!GlobalOptions || !GlobalOptions->ErrorFile)
    {
        ErrorFile = stderr;
        ErrorFile2 = NULL;
    }
    else
    {
        ErrorFile = GlobalOptions->ErrorFile;
        ErrorFile2 = stderr;
    }
    if (ErrorSeverity)
    {
        if (GlobalOptions)
        {
            GlobalOptions->ErrorCount++;
        }
        fprintf(ErrorFile, "[E%04d] %s:%d:", ErrorID, SourceFileName, SourceLine);
        if (ErrorFile2)
        {
            fprintf(ErrorFile2, "[E%04d] %s:%d:", ErrorID, SourceFileName, SourceLine);
        }
        
    }
    else
    {
        if (GlobalOptions)
        {
            GlobalOptions->WarningCount++;
        }
        fprintf(ErrorFile, "{W%04d} %s:%d:", ErrorID, SourceFileName, SourceLine);
        if (ErrorFile2)
        {
            fprintf(ErrorFile2, "{W%04d} %s:%d:", ErrorID, SourceFileName, SourceLine);
        }
    }
    ErrorMessage = "";
    for (ErrorIndex = 0; ErrorIndex < ErrorMessageCount; ErrorIndex++)
    {
        if (ErrorID == ErrorMessages[ErrorIndex].ID)
        {
            ErrorMessage = ErrorMessages[ErrorIndex].Message;
            break;
        }
    }
    switch (ArgType)
    {
    case ERROR_ARGS_S:
        fprintf(ErrorFile, ErrorMessage, StrA);
        if (ErrorFile2)
        {
            fprintf(ErrorFile2, ErrorMessage, StrA);
        }
        break;
    case ERROR_ARGS_SS:
        fprintf(ErrorFile, ErrorMessage, StrA, StrB);
        if (ErrorFile2)
        {
            fprintf(ErrorFile2, ErrorMessage, StrA, StrB);
        }
        break;
    case ERROR_ARGS_I:
        fprintf(ErrorFile, ErrorMessage, IntA);
        if (ErrorFile2)
        {
            fprintf(ErrorFile2, ErrorMessage, IntA);
        }
        break;
    case ERROR_ARGS_IS:
        fprintf(ErrorFile, ErrorMessage, IntA, StrA);
        if (ErrorFile2)
        {
            fprintf(ErrorFile2, ErrorMessage, IntA, StrA);
        }
        break;
    case ERROR_ARGS_ISS:
        fprintf(ErrorFile, ErrorMessage, IntA, StrA, StrB);
        if (ErrorFile2)
        {
            fprintf(ErrorFile2, ErrorMessage, IntA, StrA, StrB);
        }
        break;
    case ERROR_ARGS_II:
        fprintf(ErrorFile, ErrorMessage, IntA, IntB);
        if (ErrorFile2)
        {
            fprintf(ErrorFile2, ErrorMessage, IntA, IntB);
        }
        break;
    case ERROR_ARGS_IIS:
        fprintf(ErrorFile, ErrorMessage, IntA, IntB, StrA);
        if (ErrorFile2)
        {
            fprintf(ErrorFile2, ErrorMessage, IntA, IntB, StrA);
        }
        break;
    case ERROR_ARGS_III:
        fprintf(ErrorFile, ErrorMessage, IntA, IntB, IntC);
        if (ErrorFile2)
        {
            fprintf(ErrorFile2, ErrorMessage, IntA, IntB, IntC);
        }
        break;
    case ERROR_ARGS_IIII:
        fprintf(ErrorFile, ErrorMessage, IntA, IntB, IntC, IntD);
        if (ErrorFile2)
        {
            fprintf(ErrorFile2, ErrorMessage, IntA, IntB, IntC, IntD);
        }
        break;
    case ERROR_ARGS_SII:
        fprintf(ErrorFile, ErrorMessage, StrA, IntA, IntB);
        if (ErrorFile2)
        {
            fprintf(ErrorFile2, ErrorMessage, StrA, IntA, IntB);
        }
        break;
    case ERROR_ARGS_SIF:
        fprintf(ErrorFile, ErrorMessage, StrA, IntA, FloatA);
        if (ErrorFile2)
        {
            fprintf(ErrorFile2, ErrorMessage, StrA, IntA, FloatA);
        }
        break;

    case ERROR_ARGS_NONE:
    default:
        fprintf(ErrorFile, ErrorMessage);
        if (ErrorFile2)
        {
            fprintf(ErrorFile2, ErrorMessage);
        }
        break;
    }
    fprintf(ErrorFile, "\n");
    if (ErrorFile2)
    {
        fprintf(ErrorFile2, "\n");
    }

}
