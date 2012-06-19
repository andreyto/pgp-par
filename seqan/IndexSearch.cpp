// Copyright Ralf W. Grosse-Kunstleve 2002-2004. Distributed under the Boost
// Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/manage_new_object.hpp>
#include <iostream>

#include <seqan/index.h>
#include <seqan/sequence.h>

// Currently using a 5 mer since a 7 mer gives an overflow error
// during compile and a 6 mer is slower due to higher memory usage
#define MER_SIZE 5

namespace { // Avoid cluttering the global namespace.

  typedef seqan::String<seqan::AminoAcid> Peptide;
//  typedef seqan::String<char> Peptide;
//  typedef seqan::StringSet<Peptide> PepSet;
  typedef seqan::Index< Peptide,
          seqan::IndexQGram<
          seqan::UngappedShape< MER_SIZE > > > MyIndex;

  class IndexSearch
  {
    public:
      //IndexSearch(PyObject* dbStringList)
      /* Constructor takes a single string which is the Peptide database
       * to search against.
       */
      IndexSearch(const std::string& dbString)
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

          dbStr = new Peptide( dbString );
          dbIndex = new MyIndex( *dbStr );
      }
      ~IndexSearch()
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
      PyObject* find(const std::string& pepStr)
      {
          seqan::Finder<MyIndex> finder( *dbIndex );
          const Peptide pep(pepStr);
          const int pepLen = length(pep);
          PyObject* positions = PyList_New(0);

//          std::cout << "Length of db is " << seqan::length(db) << std::endl;
          while( seqan::find(finder, pep) )
          {
              int dbOffset = position(finder);
              const Peptide dbInfix = infixWithLength(*dbStr, dbOffset, pepLen);
              if (dbInfix == pepStr)
              {
                  PyObject* pos = PyInt_FromLong( dbOffset );
                  PyList_Append( positions, pos );
                  Py_DECREF( pos );
              }
          }
          return positions;
      }
    private:
//      PepSet *dbSet;
      Peptide *dbStr;
      MyIndex *dbIndex;
  };
}

BOOST_PYTHON_MODULE(IndexSearch)
{
    using namespace boost::python;
    class_<IndexSearch>("IndexSearch", init< const std::string>())
        // Add a regular member function.
        .def("find", &IndexSearch::find
            )
//                , return_value_policy<manage_new_object>() )
        ;
}

