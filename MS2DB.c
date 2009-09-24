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

// Code to parse databases in MS2DB format.

#include "CMemLeak.h"
#include "MS2DB.h"
#include "Spliced.h"
#include "expat.h"
#include "Errors.h"

#define MS2DB_BUFFER_SIZE 102400

// Macro for basic error-checking.  Assumes Cursor is set.  If the given 
// expression isn't true, we set the error flag and bail out of our current
// function.  (When Cursor->ErrorFlag is set, all our callback functions
// will terminate immediately)
#define XML_ASSERT(expr) \
    if (!(expr)) \
    {\
        Cursor->ErrorFlag = 1;\
        REPORT_ERROR_S(25, #expr);\
        return;\
    }

#define XML_ASSERT_RETVAL(expr) \
    if (!(expr)) \
    {\
        Cursor->ErrorFlag = 1;\
        REPORT_ERROR_S(25, #expr);\
        return 0;\
    }

// MS2ParseState tells us which tag we are currently inside.  
// The allowed "moves" (from tags to children) are those listed 
// in the XML schema.  However, in the interest of extensibility,
// we simply *ignore* any tags we aren't expecting to see.
typedef enum MS2ParseState
{
    evMS2DBNone = 0,
    evMS2DBDatabase,
    evMS2DBGene,
    evMS2DBGeneLocus,
    evMS2DBGeneNotes,
    evMS2DBExon,
    evMS2DBExonSequence,
    evMS2DBExonExtends,
    evMS2DBExonLinkFrom,
    evMS2DBExonMod,
    evMS2DBExonModCrossReference,
    evMS2DBGeneCrossReference
} MS2ParseState;

typedef struct MS2ParseCursor
{
    SearchInfo* Info;
    int State;
    GeneStruct* CurrentGene;
    ExonStruct* CurrentExon;
    int CurrentExonIndex;
    int CurrentLinkIndex;
    int CurrentExonSequenceIndex;
    int ErrorFlag;
    TrieNode* Root;
    XML_Parser Parser;
    int DBNumber;
} MS2ParseCursor;

// Free an MS2ParseCursor, including its attached gene.
void FreeMS2ParseCursor(MS2ParseCursor* Cursor)
{
    if (!Cursor)
    {
        return;
    }
    if (Cursor->CurrentGene)
    {
        FreeGene(Cursor->CurrentGene);
        Cursor->CurrentGene = NULL;
    }
    SafeFree(Cursor);
}

// expat callback: Handle character data in the body of a tag.
void MS2CharacterDataHandler(void* UserData, const XML_Char* String, int Length)
{
    MS2ParseCursor* Cursor;
    int NewLength;
    //
    Cursor = (MS2ParseCursor*)UserData;
    if (Cursor->ErrorFlag)
    {
        return;
    }
    switch (Cursor->State)
    {
    case evMS2DBExonSequence:
        // Incorporate this sequence into the exon sequence.
        XML_ASSERT(Cursor->CurrentExon);
        //XML_ASSERT(Cursor->CurrentExon->Sequence);
        if (!Cursor->CurrentExon->Sequence)
        {
            printf("* Warning: No sequence!?\n");
        }
        NewLength = strlen(Cursor->CurrentExon->Sequence) + Length;
        if (NewLength > Cursor->CurrentExon->Length)
        {
            REPORT_ERROR_IS(29, Cursor->CurrentExonIndex, Cursor->CurrentGene->Name);
            Cursor->ErrorFlag = 1;
            return;
        }
        strncat(Cursor->CurrentExon->Sequence, String, Length);
        break;
    default:
        break;
    }
}

// Parse attributes of a Gene tag.  
void ParseGeneAttributes(MS2ParseCursor* Cursor, const char** Attributes)
{
    int AttributeIndex;
    const char* Name;
    const char* Value;
    //
    if (Cursor->ErrorFlag)
    {
        return;
    }

    XML_ASSERT(Cursor->CurrentGene);

    for (AttributeIndex = 0; Attributes[AttributeIndex]; AttributeIndex += 2)
    {
        Name = Attributes[AttributeIndex];
        Value = Attributes[AttributeIndex + 1];
        if (!CompareStrings(Name, "Name"))
        {
            strncpy(Cursor->CurrentGene->Name, Value, GENE_NAME_LENGTH);
        }
        else if (!CompareStrings(Name, "ExonCount"))
        {
            Cursor->CurrentGene->ExonCount = atoi(Value);
        }
        else if (GlobalOptions->XMLStrictFlag)
        {
            REPORT_WARNING_ISS(28, XML_GetCurrentLineNumber(Cursor->Parser), Name, "Gene");
        }
    }
    // Allocate exons:
    XML_ASSERT(Cursor->CurrentGene->ExonCount >= 1 && Cursor->CurrentGene->ExonCount <= 10000);
    Cursor->CurrentGene->Exons = (ExonStruct*)calloc(Cursor->CurrentGene->ExonCount, sizeof(ExonStruct));
}

// Parse attributes of a Locus tag.  
void ParseLocusAttributes(MS2ParseCursor* Cursor, const char** Attributes)
{
    int AttributeIndex;
    const char* Name;
    const char* Value;
    //
    if (Cursor->ErrorFlag)
    {
        return;
    }

    XML_ASSERT(Cursor->CurrentGene);

    for (AttributeIndex = 0; Attributes[AttributeIndex]; AttributeIndex += 2)
    {
        Name = Attributes[AttributeIndex];
        Value = Attributes[AttributeIndex + 1];
        if (!CompareStrings(Name, "chromosome"))
        {
            Cursor->CurrentGene->ChromosomeNumber = atoi(Value);
        }
        else if (!CompareStrings(Name, "ForwardFlag"))
        {
            Cursor->CurrentGene->ForwardFlag = atoi(Value);
        }
        else if (GlobalOptions->XMLStrictFlag)
        {
            REPORT_WARNING_ISS(28, XML_GetCurrentLineNumber(Cursor->Parser), Name, "Locus");
        }
    }
}

// Parse attributes of an LinkFrom tag.  
// <LinkFrom Index="0" Score="3.14" AA="G" />
// If ExtendsFlag is true, this exon EXTENDS the previous one (no splicing required)
void ParseLinkFromAttributes(MS2ParseCursor* Cursor, const char** Attributes, int ExtendsFlag)
{
    int AttributeIndex;
    const char* Name;
    const char* Value;
    char EdgeAA = '\0';
    ExonEdge* Edge = NULL;
    int BackExonIndex;
    ExonStruct* Exon;
    //
    if (Cursor->ErrorFlag)
    {
        return;
    }
    Exon = Cursor->CurrentExon;
    XML_ASSERT(Exon);

    for (AttributeIndex = 0; Attributes[AttributeIndex]; AttributeIndex += 2)
    {
        Name = Attributes[AttributeIndex];
        Value = Attributes[AttributeIndex + 1];
        if (!CompareStrings(Name, "Index"))
        {
            BackExonIndex = atoi(Value);
            XML_ASSERT(BackExonIndex >= 0 && BackExonIndex < Cursor->CurrentGene->ExonCount);
            Edge = (ExonEdge*)calloc(1, sizeof(ExonEdge));
            if (ExtendsFlag)
            {
                Edge->Power = 0;
            }
            else
            {
                Edge->Power = 1;
            }
            Edge->Exon = Cursor->CurrentGene->Exons + BackExonIndex;
            Edge->Source = Exon;
            // Insert the exon into the list:
            if (Exon->BackEdgeTail)
            {
                Exon->BackEdgeTail->Next = Edge;
            }
            else
            {
                Exon->BackEdgeHead = Edge;
            }
            Exon->BackEdgeTail = Edge;
            Exon->BackEdgeCount++;
            Edge->Exon->ForwardEdgeCount++;
        }
        else if (!CompareStrings(Name, "AA"))
        {
            EdgeAA = Value[0];
        }
        else if (GlobalOptions->XMLStrictFlag)
        {
            REPORT_WARNING_ISS(28, XML_GetCurrentLineNumber(Cursor->Parser), Name, "Link");
        }
    }
    if (Edge)
    {
        Edge->AA = EdgeAA;
    }
}
// Parse attributes of an Exon tag.  
void ParseExonAttributes(MS2ParseCursor* Cursor, const char** Attributes)
{
    int AttributeIndex;
    const char* Name;
    const char* Value;
    //
    if (Cursor->ErrorFlag)
    {
        return;
    }

    XML_ASSERT(Cursor->CurrentGene);

    Cursor->CurrentExonIndex = -1; // invalidate it; the attributes will fix it
    Cursor->CurrentExon = NULL;

    // First, loop through the attributes to get the index, so we can point at the correct exon:
    for (AttributeIndex = 0; Attributes[AttributeIndex]; AttributeIndex += 2)
    {
        Name = Attributes[AttributeIndex];
        Value = Attributes[AttributeIndex + 1];
        if (!CompareStrings(Name, "index"))
        {
            Cursor->CurrentExonIndex = atoi(Value);
        }
    }
    XML_ASSERT(Cursor->CurrentExonIndex >= 0 && Cursor->CurrentExonIndex < Cursor->CurrentGene->ExonCount);

    Cursor->CurrentExon = Cursor->CurrentGene->Exons + Cursor->CurrentExonIndex;
    Cursor->CurrentExon->Gene = Cursor->CurrentGene;
    // Initialize the exon START and END to -1 (that is, not on a known chromosome):
    Cursor->CurrentExon->Start = -1;
    Cursor->CurrentExon->End = -1;
    Cursor->CurrentExonIndex = Cursor->CurrentExonIndex;

    // Now loop through and read attribute values:
    for (AttributeIndex = 0; Attributes[AttributeIndex]; AttributeIndex += 2)
    {
        Name = Attributes[AttributeIndex];
        Value = Attributes[AttributeIndex + 1];
        if (!CompareStrings(Name, "Start"))
        {
            Cursor->CurrentExon->Start = atoi(Value);
        }
        else if (!CompareStrings(Name, "End"))
        {
            Cursor->CurrentExon->End = atoi(Value);
        }
        else if (!CompareStrings(Name, "Prefix"))
        {
            strncpy(Cursor->CurrentExon->Prefix, Value, 2);
        }
        else if (!CompareStrings(Name, "Suffix"))
        {
            strncpy(Cursor->CurrentExon->Suffix, Value, 2);
        }
        else if (!CompareStrings(Name, "Index"))
        {
            ;
        }
        else if (GlobalOptions->XMLStrictFlag)
        {
            REPORT_WARNING_ISS(28, XML_GetCurrentLineNumber(Cursor->Parser), Name, "Exon");
        }
    }
    XML_ASSERT(Cursor->CurrentExonIndex >= 0 && Cursor->CurrentExonIndex < Cursor->CurrentGene->ExonCount);
}

// Parse attributes of an ExonSequence tag.  
void ParseExonSequenceAttributes(MS2ParseCursor* Cursor, const char** Attributes)
{
    int AttributeIndex;
    const char* Name;
    const char* Value;
    //
    if (Cursor->ErrorFlag)
    {
        return;
    }
    XML_ASSERT(Cursor->CurrentExon);

    for (AttributeIndex = 0; Attributes[AttributeIndex]; AttributeIndex += 2)
    {
        Name = Attributes[AttributeIndex];
        Value = Attributes[AttributeIndex + 1];
        if (!CompareStrings(Name, "Length"))
        {
            Cursor->CurrentExon->Length = atoi(Value);
        }
        else if (GlobalOptions->XMLStrictFlag)
        {
            REPORT_WARNING_ISS(28, XML_GetCurrentLineNumber(Cursor->Parser), Name, "ExonSequence");
        }
    }

    XML_ASSERT(Cursor->CurrentExon->Length >= 0 && Cursor->CurrentExon->Length < 1024*1024);
    if (Cursor->CurrentExon->Length)
    {
        Cursor->CurrentExon->Sequence = (char*)calloc(sizeof(char), Cursor->CurrentExon->Length + 1);
    }
}

// expat callback: Handle a tag and its attributes.
void MS2StartElement(void* UserData, const char* Tag, const char** Attributes)
{
    MS2ParseCursor* Cursor;
    int ExpectedTag = 0;
    //
    Cursor = (MS2ParseCursor*)UserData;
    if (Cursor->ErrorFlag)
    {
        return;
    }

    // Switch on our current state, and handle the tags we expect to see in our current state.
    // Tags we don't expect are ignored (i.e. new tags can be added without breaking the parser)
    switch (Cursor->State)
    {
    case evMS2DBNone:
        if (!CompareStrings(Tag, "Database"))
        {
            ExpectedTag = 1;
            Cursor->State = evMS2DBDatabase;
            // ignore database attributes for now
        }
        break;
    case evMS2DBDatabase:
        if (!CompareStrings(Tag, "Gene"))
        {
            ExpectedTag = 1;
            XML_ASSERT(!Cursor->CurrentGene);
            Cursor->State = evMS2DBGene;
            Cursor->CurrentGene = (GeneStruct*)calloc(1, sizeof(GeneStruct));
            Cursor->CurrentExonIndex = 0;
            Cursor->CurrentGene->ChromosomeNumber = -1;
            Cursor->CurrentGene->ForwardFlag = 1; // default
            ParseGeneAttributes(Cursor, Attributes);
        }
        if (!CompareStrings(Tag, "Locus"))
        {
            ExpectedTag = 1;
            Cursor->State = evMS2DBGeneLocus;
            ParseLocusAttributes(Cursor, Attributes);
        }
        break;
    case evMS2DBGene:
        XML_ASSERT(Cursor->CurrentGene);
        if (!CompareStrings(Tag, "Exon"))
        {
            ExpectedTag = 1;
            Cursor->State = evMS2DBExon;
            ParseExonAttributes(Cursor, Attributes);
        }
        if (!CompareStrings(Tag, "CrossReference"))
        {
            // We don't do anything with the attributes, but cross-references
            // are "expected", so we don't raise a warning.
            ExpectedTag = 1;
            Cursor->State = evMS2DBGeneCrossReference;
        }
        break;
    case evMS2DBGeneCrossReference:
        if (!CompareStrings(Tag, "CRExons"))
        {
            ExpectedTag = 1;
        }
        break;
    case evMS2DBExon:
        XML_ASSERT(Cursor->CurrentExon);
        if (!CompareStrings(Tag, "ExonSequence"))
        {
            ExpectedTag = 1;
            Cursor->State = evMS2DBExonSequence;
            ParseExonSequenceAttributes(Cursor, Attributes);
        }
        if (!CompareStrings(Tag, "ExtendsExon"))
        {
            ExpectedTag = 1;
            // Don't change states, ExtendsExon has no body
            ParseLinkFromAttributes(Cursor, Attributes, 1);
        }
        if (!CompareStrings(Tag, "LinkFrom"))
        {
            ExpectedTag = 1;
            // Don't change states, LinkFrom has no body
            ParseLinkFromAttributes(Cursor, Attributes, 0);
        }
        break;
    default:
        break;
    }
    if (!ExpectedTag)
    {
        REPORT_ERROR_IS(27, XML_GetCurrentLineNumber(Cursor->Parser), Tag);
    }
}

// Confirm that this gene is, indeed, searchable.
int IntegrityCheckXMLGene(MS2ParseCursor* Cursor)
{
    int ExonIndex;
    int EdgeIndex;
    ExonEdge* Edge;
    ExonEdge* PrevEdge;
    GeneStruct* Gene;
    ExonStruct* Exon;
    //
    Gene = Cursor->CurrentGene;
    XML_ASSERT_RETVAL(Gene);
    for (ExonIndex = 0; ExonIndex < Cursor->CurrentGene->ExonCount; ExonIndex++)
    {
        // Confirm that we did, in fact, observe this exon:
        Exon = Gene->Exons + ExonIndex;
        if (!Exon->Gene)
        {
            printf("* Error: Exon '%d' from Gene '%s' not present!\n", ExonIndex, Gene->Name);
            return 0;
        }
    }
    // All exons have been initialized.  Now let's fix up the backward edges from each exon.
    // We MOVE the backward edges from the linked list Exon->BackEdgeHead->Next->...->Exon->BackEdgeTail
    // into an array, Exon->BackwardEdges.
    for (ExonIndex = 0; ExonIndex < Cursor->CurrentGene->ExonCount; ExonIndex++)
    {
        Exon = Gene->Exons + ExonIndex;
        PrevEdge = NULL;
        if (Exon->BackEdgeCount)
        {
            Exon->BackwardEdges = (ExonEdge*)calloc(Exon->BackEdgeCount, sizeof(ExonEdge));
            for (EdgeIndex = 0, Edge = Exon->BackEdgeHead; Edge; EdgeIndex++, Edge = Edge->Next)
            {
                memcpy(Exon->BackwardEdges + EdgeIndex, Edge, sizeof(ExonEdge));
                SafeFree(PrevEdge);
                PrevEdge = Edge;
            }
            SafeFree(PrevEdge);
            Exon->BackEdgeHead = NULL;
            Exon->BackEdgeTail = NULL;
        }
        // Allocate forward-edge array:
        Exon->ForwardEdges = (ExonEdge*)calloc(Exon->ForwardEdgeCount, sizeof(ExonEdge));
    }
    // Finally, we'll set the forward edges from each exon:
    SetExonForwardEdges(Cursor->CurrentGene);
    return 1;
}

void MS2EndElement(void* UserData, const char* Tag)
{
    MS2ParseCursor* Cursor;
    int Result;
    //
    Cursor = (MS2ParseCursor*)UserData;
    if (Cursor->ErrorFlag)
    {
        return;
    }
    //printf("End tag '%s', current state %d\n", Tag, Cursor->State);
    switch (Cursor->State)
    {
    case evMS2DBDatabase:
        if (!CompareStrings(Tag, "Database"))
        {
            Cursor->State = evMS2DBNone;
        }
        break;
    case evMS2DBGene:
        if (!CompareStrings(Tag, "Gene"))
        {
            Cursor->State = evMS2DBDatabase;
            // We search a gene immediately after we finish parsing it.  (Note that by the 
            // time control returns to SearchMS2DB, we may have shot through 10 genes!)
            Result = IntegrityCheckXMLGene(Cursor);
            if (Result)
            {
                SearchSplicableGene(Cursor->Info, Cursor->CurrentGene);
            }
            FreeGene(Cursor->CurrentGene);
            Cursor->CurrentGene = NULL;
            Cursor->Info->RecordNumber++;
        }
        break;
    case evMS2DBExon:
        if (!CompareStrings(Tag, "Exon"))
        {
            Cursor->State = evMS2DBGene;
            //Cursor->CurrentExonIndex++;
        }
        break;
    case evMS2DBExonSequence:
        if (!CompareStrings(Tag, "ExonSequence"))
        {
            Cursor->State = evMS2DBExon;
        }
        break;
    case evMS2DBGeneLocus:
        if (!CompareStrings(Tag, "Locus"))
        {
            Cursor->State = evMS2DBGene;
        }
        break;
    case evMS2DBGeneCrossReference:
        if (!CompareStrings(Tag, "CrossReference"))
        {
            Cursor->State = evMS2DBGene;
        }
        break;
    default:
        printf("* Error: End-tag '%s' not handled from state %d\n", Tag, Cursor->State);
        Cursor->ErrorFlag = 1;
        break;
    }
}

void SearchMS2DB(SearchInfo* Info)
{
    FILE* DBFile;
    XML_Parser Parser = NULL;
    int ParseUserData = 0;
    int XMLParseResult;
    int BytesRead;
    int DoneFlag = 0;
    int FilePos = 0;
    void* XMLBuffer;
    MS2ParseCursor* Cursor;
    int Error;
    //
    DBFile = Info->DB->DBFile;
    if (!DBFile)
    {
        printf("** Error: Unable to open database file '%s'\n", Info->DB->FileName);
        return;
    }
    fseek(DBFile, 0, 0);
    AllocSpliceStructures();
    Cursor = (MS2ParseCursor*)calloc(sizeof(MS2ParseCursor), 1);
    Cursor->Info = Info;
    Parser = XML_ParserCreate(NULL);
    Cursor->Parser = Parser;
    XML_SetUserData(Parser, Cursor);
    XML_SetElementHandler(Parser, MS2StartElement, MS2EndElement);
    XML_SetCharacterDataHandler(Parser, MS2CharacterDataHandler);

    while (!DoneFlag)
    {
        // Get a buffer (parser handles the memory):
        XMLBuffer = XML_GetBuffer(Parser, sizeof(char) * MS2DB_BUFFER_SIZE);
        if (!XMLBuffer)
        {
            printf("* Error: Unable to get XML buffer of size %d\n", MS2DB_BUFFER_SIZE);
            break;
        }

        // Read into the buffer:
        BytesRead = ReadBinary(XMLBuffer, sizeof(char), MS2DB_BUFFER_SIZE, DBFile);
        if (!BytesRead)
        {
            // We'll call XML_Parse once more, this time with DoneFlag set to 1. 
            DoneFlag = 1;
        }

        // Parse this block o' text:
        XMLParseResult = XML_Parse(Parser, XMLBuffer, BytesRead, DoneFlag);
        if (!XMLParseResult)
        {
            printf("XML Parse error - file position ~%d\n", FilePos);
            Error = XML_GetErrorCode(Parser);
            printf("Error code %d description '%s'\n", Error, XML_ErrorString(Error));
        }

        // If Cursor->ErrorFlag is set, then the file isn't valid!  Error out
        // now, since recovery could be difficult.
        if (Cursor->ErrorFlag)
        {
            break;
        }
        FilePos += BytesRead;
    }

    XML_ParserFree(Parser);
    SafeFree(Cursor);
}


