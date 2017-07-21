/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef UNITTESTFIELDUTILS_H
#define UNITTESTFIELDUTILS_H

#include <gtest/gtest.h>
#include "UnitTestUtils.h"
#include "UnitTestKokkosUtils.h"

namespace unit_test_utils {

double field_norm(const ScalarFieldType& field, const stk::mesh::BulkData& bulk, stk::mesh::Selector selector);

}

#endif /* UNITTESTFIELDUTILS_H */
