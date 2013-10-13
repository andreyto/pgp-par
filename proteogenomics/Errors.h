#ifndef ERRORS_H
#define ERRORS_H
//
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

void AssertionFailed(char* Assertion, char* FileName, int LineNumber);
void InitErrors();
// ReportError is not to be called directly.  Use the REPORT_ERROR and REPORT_WARNING macros.
void ReportError(int ErrorSeverity, int ErrorID, int SourceLine, char* SourceFileName, int Args,
                 const char* StrA, const char* StrB, const char* StrC, 
                 int IntA, int IntB, int IntC, int IntD,
                 float FloatA, float FloatB);

#define ERROR_ARGS_NONE 0
#define ERROR_ARGS_S 1
#define ERROR_ARGS_I 2
#define ERROR_ARGS_II 3
#define ERROR_ARGS_III 4
#define ERROR_ARGS_IIII 5
#define ERROR_ARGS_IS 6
#define ERROR_ARGS_IIS 7
#define ERROR_ARGS_IIIS 8
#define ERROR_ARGS_SS 9
#define ERROR_ARGS_ISS 10
#define ERROR_ARGS_SII 11
#define ERROR_ARGS_SIF 12

#define REPORT_ERROR(ErrorID) ReportError(1, ErrorID, __LINE__, __FILE__, ERROR_ARGS_NONE, NULL, NULL, NULL, 0, 0, 0, 0, 0.0, 0.0);
#define REPORT_ERROR_S(ErrorID, StrA) ReportError(1, ErrorID, __LINE__, __FILE__, ERROR_ARGS_S, StrA, NULL, NULL, 0, 0, 0, 0, 0.0, 0.0);
#define REPORT_ERROR_SS(ErrorID, StrA, StrB) ReportError(1, ErrorID, __LINE__, __FILE__, ERROR_ARGS_SS, StrA, StrB, NULL, 0, 0, 0, 0, 0.0, 0.0);
#define REPORT_ERROR_I(ErrorID, IntA) ReportError(1, ErrorID, __LINE__, __FILE__, ERROR_ARGS_I, NULL, NULL, NULL, IntA, 0, 0, 0, 0.0, 0.0);
#define REPORT_ERROR_II(ErrorID, IntA, IntB) ReportError(1, ErrorID, __LINE__, __FILE__, ERROR_ARGS_II, NULL, NULL, NULL, IntA, IntB, 0, 0, 0.0, 0.0);
#define REPORT_ERROR_III(ErrorID, IntA, IntB, IntC) ReportError(1, ErrorID, __LINE__, __FILE__, ERROR_ARGS_III, NULL, NULL, NULL, IntA, IntB, IntC, 0, 0.0, 0.0);
#define REPORT_ERROR_IIII(ErrorID, IntA, IntB, IntC, IntD) ReportError(1, ErrorID, __LINE__, __FILE__, ERROR_ARGS_IIII, NULL, NULL, NULL, IntA, IntB, IntC, IntD, 0.0, 0.0);
#define REPORT_ERROR_IS(ErrorID, IntA, StrA) ReportError(1, ErrorID, __LINE__, __FILE__, ERROR_ARGS_IS, StrA, NULL, NULL, IntA, 0, 0, 0, 0.0, 0.0);
#define REPORT_ERROR_IIS(ErrorID, IntA, IntB, StrA) ReportError(1, ErrorID, __LINE__, __FILE__, ERROR_ARGS_IIS, StrA, NULL, NULL, IntA, IntB, 0, 0, 0.0, 0.0);
#define REPORT_ERROR_IIIS(ErrorID, IntA, IntB, IntC, StrA) ReportError(1, ErrorID, __LINE__, __FILE__, ERROR_ARGS_IIIS, StrA, NULL, NULL, IntA, IntB, IntC, 0, 0.0, 0.0);
#define REPORT_ERROR_ISS(ErrorID, IntA, StrA, StrB) ReportError(1, ErrorID, __LINE__, __FILE__, ERROR_ARGS_ISS, StrA, StrB, NULL, IntA, 0, 0, 0, 0.0, 0.0);
#define REPORT_ERROR_SII(ErrorID, StrA, IntA, IntB) ReportError(1, ErrorID, __LINE__, __FILE__, ERROR_ARGS_SII, StrA, NULL, NULL, IntA, IntB, 0, 0, 0.0, 0.0);
#define REPORT_ERROR_SIF(ErrorID, StrA, IntA, FloatA) ReportError(1, ErrorID, __LINE__, __FILE__, ERROR_ARGS_SII, StrA, NULL, NULL, IntA, 0, 0, 0, FloatA, 0.0);

#define REPORT_WARNING(ErrorID) ReportError(0, ErrorID, __LINE__, __FILE__, ERROR_ARGS_NONE, NULL, NULL, NULL, 0, 0, 0, 0, 0.0, 0.0);
#define REPORT_WARNING_S(ErrorID, StrA) ReportError(0, ErrorID, __LINE__, __FILE__, ERROR_ARGS_S, StrA, NULL, NULL, 0, 0, 0, 0, 0.0, 0.0);
#define REPORT_WARNING_SS(ErrorID, StrA, StrB) ReportError(0, ErrorID, __LINE__, __FILE__, ERROR_ARGS_SS, StrA, StrB, NULL, 0, 0, 0, 0, 0.0, 0.0);
#define REPORT_WARNING_I(ErrorID, IntA) ReportError(0, ErrorID, __LINE__, __FILE__, ERROR_ARGS_I, NULL, NULL, NULL, IntA, 0, 0, 0, 0.0, 0.0);
#define REPORT_WARNING_II(ErrorID, IntA, IntB) ReportError(0, ErrorID, __LINE__, __FILE__, ERROR_ARGS_II, NULL, NULL, NULL, IntA, IntB, 0, 0, 0.0, 0.0);
#define REPORT_WARNING_III(ErrorID, IntA, IntB, IntC) ReportError(1, ErrorID, __LINE__, __FILE__, ERROR_ARGS_III, NULL, NULL, NULL, IntA, IntB, IntC, 0, 0.0, 0.0);
#define REPORT_WARNING_IIII(ErrorID, IntA, IntB, IntC, IntD) ReportError(1, ErrorID, __LINE__, __FILE__, ERROR_ARGS_IIII, NULL, NULL, NULL, IntA, IntB, IntC, IntD, 0.0, 0.0);
#define REPORT_WARNING_IS(ErrorID, IntA, StrA) ReportError(0, ErrorID, __LINE__, __FILE__, ERROR_ARGS_IS, StrA, NULL, NULL, IntA, 0, 0, 0, 0.0, 0.0);
#define REPORT_WARNING_IIS(ErrorID, IntA, IntB, StrA) ReportError(0, ErrorID, __LINE__, __FILE__, ERROR_ARGS_IIS, StrA, NULL, NULL, IntA, IntB, 0, 0, 0.0, 0.0);
#define REPORT_WARNING_IIIS(ErrorID, IntA, IntB, IntC, StrA) ReportError(0, ErrorID, __LINE__, __FILE__, ERROR_ARGS_IIIS, StrA, NULL, NULL, IntA, IntB, IntC, 0, 0.0, 0.0);
#define REPORT_WARNING_ISS(ErrorID, IntA, StrA, StrB) ReportError(0, ErrorID, __LINE__, __FILE__, ERROR_ARGS_ISS, StrA, StrB, NULL, IntA, 0, 0, 0, 0.0, 0.0);
#define REPORT_WARNING_SII(ErrorID, StrA, IntA, IntB) ReportError(0, ErrorID, __LINE__, __FILE__, ERROR_ARGS_SII, StrA, NULL, NULL, IntA, IntB, 0, 0, 0.0, 0.0);
#define REPORT_WARNING_SIF(ErrorID, StrA, IntA, FloatA) ReportError(0, ErrorID, __LINE__, __FILE__, ERROR_ARGS_SII, StrA, NULL, NULL, IntA, 0, 0, 0, FloatA, 0.0);

#endif // ERRORS_H
