// Copyright Ralf W. Grosse-Kunstleve 2002-2004. Distributed under the Boost
// Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/manage_new_object.hpp>
#include <iostream>

#include <seqan/find.h>
#include <seqan/sequence.h>

namespace { // Avoid cluttering the global namespace.

  //typedef seqan::String<seqan::AminoAcid> Peptide;
  typedef seqan::String<char> Peptide;

  class ExactSearch
  {
    public:
      ExactSearch(){}

      PyObject* find(const std::string& dbStr, const std::string& pepStr)
      {
          Peptide db(dbStr);
          seqan::Finder<Peptide> finder(db);
          const Peptide pep(pepStr);
          seqan::Pattern<Peptide, seqan::Horspool> pattern(pep);

          PyObject* positions = PyList_New(0);

//          std::cout << "Length of db is " << seqan::length(db) << std::endl;
          while( seqan::find(finder, pattern) )
          {
              PyObject* pos = PyInt_FromLong( position( finder ));
              PyList_Append( positions, pos );
              Py_DECREF( pos );
          }
          return positions;
      }
    private:
  };
}

BOOST_PYTHON_MODULE(ExactSearch)
{
    using namespace boost::python;
    class_<ExactSearch>("ExactSearch")
        // Add a regular member function.
        .def("find", &ExactSearch::find
            )
//                , return_value_policy<manage_new_object>() )
        ;
}

