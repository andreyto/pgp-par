#include <Python.h>

//#include <iostream>

#include <seqan/index.h>
#include <seqan/sequence.h>

// Currently using a 5 mer since a 7 mer gives an overflow error
// during compile and a 6 mer is slower due to higher memory usage
#define MER_SIZE 5

namespace
{				// Avoid cluttering the global namespace.

  typedef seqan::String < seqan::AminoAcid > Peptide;
//  typedef seqan::String<char> Peptide;
//  typedef seqan::StringSet<Peptide> PepSet;
  typedef seqan::Index < Peptide,
    seqan::IndexQGram < seqan::UngappedShape < MER_SIZE > > >MyIndex;

  class IndexSearch
  {
  public:
    //IndexSearch(PyObject* dbStringList)
    /* Constructor takes a single string which is the Peptide database
     * to search against.
     */
    IndexSearch (const std::string & dbString)
    {
      /*
         PyObject *iter = PyObject_GetIter( dbStringList );
         PyObject *dbSeq;
         dbSet = new PepSet();
         resize( dbSet, PyList_Size( dbStringList ) );

         int i = 0;
         while ( dbSeq = PyIter_Next( iter ) )
         {
         Peptide nextPep( PyString_AsString( dbSeq ));
         dbSet[i++] = nextPep;
         }
       */

      dbStr = new Peptide (dbString);
      dbIndex = new MyIndex (*dbStr);
    }
     ~IndexSearch ()
    {
      delete dbIndex;
      delete dbStr;
    }

    /* find takes a single peptide string which is searched against the db.
     * It locates the 1st kmer in the database using seqan::find. Finally,
     * it does an exact string compare of the search string on all the kmer
     * positions in the database to find the full matches, which are returned
     * in a python list.
     */
    PyObject *find (const std::string & pepStr)
    {
      seqan::Finder < MyIndex > finder (*dbIndex);
      const Peptide pep (pepStr);
      const int pepLen = length (pep);
      PyObject *positions = PyList_New (0);

//          std::cout << "Length of db is " << seqan::length(db) << std::endl;
      while (seqan::find (finder, pep))
	{
	  int dbOffset = position (finder);
	  const Peptide dbInfix = infixWithLength (*dbStr, dbOffset, pepLen);
	  if (dbInfix == pepStr)
	    {
	      PyObject *pos = PyInt_FromLong (dbOffset);
	      PyList_Append (positions, pos);
	      Py_DECREF (pos);
	    }
	}
      return positions;
    }
  private:
//      PepSet *dbSet;
    Peptide * dbStr;
    MyIndex *dbIndex;
  };
}

#ifdef __cplusplus
extern "C" {
#endif

const char *caps_name = "IndexSearch";

static void del_IndexSearch(PyObject *pcaps)
{
	IndexSearch * oldind = static_cast<IndexSearch *>(PyCapsule_GetPointer(pcaps, caps_name));
    delete oldind;
}

static PyObject *new_IndexSearch(PyObject *, PyObject* args)
{
    char *dbString = 0;
    int ok = PyArg_ParseTuple(args,"s",&dbString);
    if(!ok) return NULL;

    try {
	    //Capsule is the standard way to pass opaque C++ pointers
        //around starting with Python 2.7
        return PyCapsule_New(new IndexSearch(dbString), caps_name, del_IndexSearch);
    }
    catch (...) {

        PyErr_SetString( PyExc_ValueError, 
                     "Unable to create IndexSearch");
        return NULL;    // trigger exception
    }
}



static PyObject *IndexSearch_find(PyObject *, PyObject* args)
{
    // First, get the PyCObject from the args tuple
    PyObject *pcaps = 0;
    char *pepStr = 0;
    int ok = PyArg_ParseTuple( args, "Os", &pcaps,&pepStr);
    //"O" is for Object
    if(!ok) return NULL;

	IndexSearch * thisind = static_cast<IndexSearch *>(PyCapsule_GetPointer(pcaps, caps_name));

    return thisind->find(pepStr);
}

static PyMethodDef indexSearchMethods[] = 
{
    { "new_IndexSearch", new_IndexSearch, 
        METH_VARARGS, 
      "new_IndexSearch(str)->new IndexSearch object"},
    { "IndexSearch_find", IndexSearch_find, 
      METH_VARARGS, 
      "IndexSearch_find(IndexSearch,str) -> list"},

    {NULL,NULL,0,NULL}
};

#ifdef __cplusplus
}
#endif

PyMODINIT_FUNC
init_IndexSearch()
{
       Py_InitModule3(
            "_IndexSearch",   // name of the module
            indexSearchMethods,  // name of the method table
            "C++ IndexSearch class"); // doc string for module
}

