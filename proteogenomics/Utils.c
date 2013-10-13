
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
#include "Utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <stdarg.h>

// From high to low
int CompareFloats(const float* a, const float* b)
{
    if (*a > *b)
    {
        return -1;
    }
    if (*a < *b)
    {
        return 1;
    }
    return 0;
}

// From high to low
int CompareInts(const int* a, const int* b)
{
    if (*a > *b)
    {
        return -1;
    }
    if (*a < *b)
    {
        return 1;
    }
    return 0;
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

void ParseFileByLines(FILE* File, FileLineParser LineParser, void* ParseData, int ProcessCommentLines)
{
    char LineBuffer[MAX_LINE_LENGTH];
    char TextBuffer[BUFFER_SIZE * 2];
    int LineNumber = 0;
    int FilePos;
    int NewFilePos = 0;
    int BytesToRead;
    int BufferEnd = 0;
    int BytesRead;
    int BufferPos = 0;
    int KeepParsingFlag;
    //
    if (!File)
    {
        return;
    }
    NewFilePos = ftell(File);
    while (1)
    {
        FilePos = NewFilePos;
        BytesToRead = BUFFER_SIZE - BufferEnd;
        BytesRead = ReadBinary(TextBuffer + BufferEnd, sizeof(char), BytesToRead, File);
        BufferEnd += BytesRead;
        TextBuffer[BufferEnd] = '\0';
        if (BufferPos == BufferEnd)
        { 
            // We're done!
            break;
        }
        BufferPos = CopyBufferLine(TextBuffer, BufferPos, BufferEnd, LineBuffer, 0);
        if (!BufferPos)
        {
            // We encountered a null character.  Force advance:
            BufferPos++;
        }
        NewFilePos = FilePos + BufferPos;
        
        LineNumber += 1;
        // Now, move the remaining text to the start of the buffer:
        memmove(TextBuffer, TextBuffer + BufferPos, BufferEnd - BufferPos);
        BufferEnd -= BufferPos;
        BufferPos = 0;
        // Now, process this line of text!
        // Skip empty lines:
        if (!LineBuffer[0])
        {
            continue;
        }
        // Skip comment lines:
        if (LineBuffer[0] == '#' && !ProcessCommentLines)
        {
            continue;
        }
        KeepParsingFlag = LineParser(LineNumber, FilePos, LineBuffer, ParseData);
        if (!KeepParsingFlag)
        {
            break;
        }
    }
}

#define FORCE_UPPER(X) X = ((X) >= 'A' && (X) <= 'Z' ? (X) + 'a' - 'A' : (X));

// Case-insensitive string comparison.  Returns -1 if A<B, 1 if A>B, 0 if same.
int CompareStrings(const char* StringA, const char* StringB)
{
    const char* CharA;
    const char* CharB;
    char A;
    char B;
    CharA = StringA;
    CharB = StringB;
    while (1)
    {
        if (!*CharA && !*CharB)
        {
            return 0;
        }
        A = *CharA;
        B = *CharB;
        FORCE_UPPER(A);
        FORCE_UPPER(B);
        //if (isupper(A)) A = ConvertToLower(A);
        //if (isupper(B)) B = ConvertToLower(B);
        if (A < B)
        {
            return 1;
        }
        if (A > B)
        {
            return -1;
        }
        CharA++;
        CharB++;
    }
}

#ifdef __ppc__
// Reads a little endian binary file for a big endian system
size_t ReadBinary(void* Buffer, size_t ItemSize, size_t ItemCount, FILE* File) 
{
    size_t ItemIndex;
    size_t ByteIndex;
    unsigned char SwapValue;
    char* CharBuffer;
    int BytesRead;

    BytesRead = fread(Buffer, size, MemberCount, File); // raw fread

    CharBuffer = (char*)Buffer;

    for (ItemIndex = 0; ItemIndex < ItemCount; ItemIndex++) 
    {
        for (ByteIndex = 0; ByteIndex < ItemSize >> 1; ByteIndex++)
        {
            // Swap the first and last bytes, then bytes 1 and max - 1, etc.
            SwapValue = CharBuffer[size * ItemSize + ByteIndex];
            CharBuffer[size * ItemIndex + ByteIndex] = CharBuffer[ItemSize * ItemIndex + ItemSize - ByteIndex - 1];
            CharBuffer[size * ItemIndex + ItemSize - ByteIndex - 1] = SwapValue;
        }
    }
    return BytesRead;
}

// We're on a big-endian system, and we must write out a little-endian file.
size_t WriteBinary(void* Buffer, size_t ItemSize, size_t ItemCount, FILE* File)
{
    char ItemBuffer[256];
    int ItemIndex;
    int ByteIndex;
    char* CharBuffer = (char*)Buffer;
    //
    // Write a byte-swapped version of each item to ItemBuffer, then write ItemBuffer to disk.
    for (ItemIndex = 0; ItemIndex < ItemCount; ItemIndex++) 
    {
        for (ByteIndex = 0; ByteIndex < ItemSize; ByteIndex++)
        {
            ItemBuffer[ItemSize - ByteIndex - 1] = CharBuffer[ItemSize * ItemIndex + ByteIndex]
            fwrite(ItemBuffer, ItemSize, 1, File); // raw fwrite
        }
    }
}

#else
#define ReadBinary fread
#define WriteBinary fwrite
#endif

char TranslateCodon(char* DNA)
{
    switch (DNA[0])
    {
    case 'T':
    case 't':
        switch (DNA[1])
        {
        case 'T':
        case 't':
            switch (DNA[2])
            {
            case 'T':
            case 't':
                return 'F';
            case 'C':
            case 'c':
                return 'F';
            case 'A':
            case 'a':
                return 'L';
            case 'G':
            case 'g':
                return 'L';
            }
            break;
        case 'C':
        case 'c':
            switch (DNA[2])
            {
            case 'T':
            case 't':
                return 'S';
            case 'C':
            case 'c':
                return 'S';
            case 'A':
            case 'a':
                return 'S';
            case 'G':
            case 'g':
                return 'S';
            }
            break;
        case 'A':
        case 'a':
            switch (DNA[2])
            {
            case 'T':
            case 't':
                return 'Y';
            case 'C':
            case 'c':
                return 'Y';
            case 'A':
            case 'a':
                return 'X';
            case 'G':
            case 'g':
                return 'X';
            }
            break;
        case 'G':
        case 'g':
            switch (DNA[2])
            {
            case 'T':
            case 't':
                return 'C';
            case 'C':
            case 'c':
                return 'C';
            case 'A':
            case 'a':
                return 'X';
            case 'G':
            case 'g':
                return 'W';
            }
            break;
        }
        break;
    case 'C':
    case 'c':
        switch (DNA[1])
        {
        case 'T':
        case 't':
            switch (DNA[2])
            {
            case 'T':
            case 't':
                return 'L';
            case 'C':
            case 'c':
                return 'L';
            case 'A':
            case 'a':
                return 'L';
            case 'G':
            case 'g':
                return 'L';
            }
            break;
        case 'C':
        case 'c':
            switch (DNA[2])
            {
            case 'T':
            case 't':
                return 'P';
            case 'C':
            case 'c':
                return 'P';
            case 'A':
            case 'a':
                return 'P';
            case 'G':
            case 'g':
                return 'P';
            }
            break;
        case 'A':
        case 'a':
            switch (DNA[2])
            {
            case 'T':
            case 't':
                return 'H';
            case 'C':
            case 'c':
                return 'H';
            case 'A':
            case 'a':
                return 'Q';
            case 'G':
            case 'g':
                return 'Q';
            }
            break;
        case 'G':
        case 'g':
            switch (DNA[2])
            {
            case 'T':
            case 't':
                return 'R';
            case 'C':
            case 'c':
                return 'R';
            case 'A':
            case 'a':
                return 'R';
            case 'G':
            case 'g':
                return 'R';
            }
            break;
        }
        break;
    case 'A':
    case 'a':
        switch (DNA[1])
        {
        case 'T':
        case 't':
            switch (DNA[2])
            {
            case 'T':
            case 't':
                return 'I';
            case 'C':
            case 'c':
                return 'I';
            case 'A':
            case 'a':
                return 'I';
            case 'G':
            case 'g':
                return 'M';
            }
            break;
        case 'C':
        case 'c':
            switch (DNA[2])
            {
            case 'T':
            case 't':
                return 'T';
            case 'C':
            case 'c':
                return 'T';
            case 'A':
            case 'a':
                return 'T';
            case 'G':
            case 'g':
                return 'T';
            }
            break;
        case 'A':
        case 'a':
            switch (DNA[2])
            {
            case 'T':
            case 't':
                return 'N';
            case 'C':
            case 'c':
                return 'N';
            case 'A':
            case 'a':
                return 'K';
            case 'G':
            case 'g':
                return 'K';
            }
            break;
        case 'G':
        case 'g':
            switch (DNA[2])
            {
            case 'T':
            case 't':
                return 'S';
            case 'C':
            case 'c':
                return 'S';
            case 'A':
            case 'a':
                return 'R';
            case 'G':
            case 'g':
                return 'R';
            }
            break;
        }
        break;
    case 'G':
    case 'g':
        switch (DNA[1])
        {
        case 'T':
        case 't':
            switch (DNA[2])
            {
            case 'T':
            case 't':
                return 'V';
            case 'C':
            case 'c':
                return 'V';
            case 'A':
            case 'a':
                return 'V';
            case 'G':
            case 'g':
                return 'V';
            }
            break;
        case 'C':
        case 'c':
            switch (DNA[2])
            {
            case 'T':
            case 't':
                return 'A';
            case 'C':
            case 'c':
                return 'A';
            case 'A':
            case 'a':
                return 'A';
            case 'G':
            case 'g':
                return 'A';
            }
            break;
        case 'A':
        case 'a':
            switch (DNA[2])
            {
            case 'T':
            case 't':
                return 'D';
            case 'C':
            case 'c':
                return 'D';
            case 'A':
            case 'a':
                return 'E';
            case 'G':
            case 'g':
                return 'E';
            }
            break;
        case 'G':
        case 'g':
            switch (DNA[2])
            {
            case 'T':
            case 't':
                return 'G';
            case 'C':
            case 'c':
                return 'G';
            case 'A':
            case 'a':
                return 'G';
            case 'G':
            case 'g':
                return 'G';
            }
            break;
        }
        break;
    }
    return 'X';
}

void WriteReverseComplement(char* Source, char* Destination)
{
    char* A;
    char* B;
    A = Source;
    while (*A)
    {
        A++;
    }
    A--;
    B = Destination;
    while (A >= Source)
    {
        switch (*A)
        {
        case 'C':
        case 'c':
            *B = 'G';
            break;
        case 'G':
        case 'g':
            *B = 'C';
            break;
        case 'A':
        case 'a':
            *B = 'T';
            break;
        case 'T':
        case 't':
            *B = 'A';
            break;
        }
        A--;
        B++;
    }
}

// Reverse a null-terminated string in place:
void ReverseString(char* String)
{
    char* A;
    char* Z;
    char Temp;
    int Len;
    if (!String)
    {
        return;
    }
    Len = strlen(String);
    if (!Len)
    {
        return;
    }
    Z = String + Len - 1;
    A = String;
    while (A < Z)
    {
        Temp = *Z;
        *Z = *A;
        *A = Temp;
        A++;
        Z--;
    }
}

float GetMedian(float* Values, int ValueCount)
{
    qsort(Values, ValueCount, sizeof(float), (QSortCompare)CompareFloats);
    if (ValueCount % 2)
    {
        return Values[ValueCount / 2];
    }
    else
    {
        return (Values[ValueCount / 2] + Values[(ValueCount / 2) - 1]) / (float)2.0;
    }
}
