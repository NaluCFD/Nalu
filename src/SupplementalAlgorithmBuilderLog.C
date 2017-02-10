/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#include <SupplementalAlgorithmBuilderLog.h>
#include <NaluEnv.h>

#include <map>
#include <set>
#include <string>
#include <vector>


namespace sierra{
namespace nalu{
//--------------------------------------------------------------------------
SuppAlgBuilderLog&
SuppAlgBuilderLog::self()
{
  static SuppAlgBuilderLog instance;
  return instance;
}
//--------------------------------------------------------------------------
void
SuppAlgBuilderLog::add_valid_name(std::string suppAlgTypeName, std::string name)
{
  validSuppAlgNames_[suppAlgTypeName].insert(name);
}
//--------------------------------------------------------------------------
void
SuppAlgBuilderLog::add_built_name(std::string suppAlgTypeName, std::string name)
{
  builtSuppAlgNames_[suppAlgTypeName].insert(name);
}
//--------------------------------------------------------------------------
bool
SuppAlgBuilderLog::print_invalid_supp_alg_names(
  std::string suppAlgTypeName,
  const std::map<std::string, std::vector<std::string>>& srcTermsMap )
{
  auto it = srcTermsMap.find(suppAlgTypeName);
  if (it == srcTermsMap.end()) {
    return false;
  }
  const std::vector<std::string>& inputFileNames = it->second;
  const std::set<std::string> implementedSuppAlgNames = valid_supp_alg_names(suppAlgTypeName);

  std::vector<std::string> badNames;
  for (const auto& inputFileName : inputFileNames) {
    if (implementedSuppAlgNames.find(inputFileName) == implementedSuppAlgNames.end()) {
      badNames.push_back(inputFileName);
    }
  }

  for (const auto& name : badNames) {
    NaluEnv::self().naluOutputP0() << "Error: No Supplemental Algorithm with name `"
                                   << name
                                   << "' implemented"
                                   << std::endl;
  }

  bool isOK = badNames.empty();

  return isOK;
}
//--------------------------------------------------------------------------
void
SuppAlgBuilderLog::print_valid_supp_alg_names(std::string suppAlgTypeName )
{
  const auto implementedSuppAlgNames = valid_supp_alg_names(suppAlgTypeName);

  std::string msgList = "";
  for (const auto& name : implementedSuppAlgNames) {
    msgList += "`" + name + "', ";
  }
  msgList = msgList.substr(0, msgList.size()-2);

  NaluEnv::self().naluOutputP0()
      << "Valid Supplemental Algorithm names for "
      << suppAlgTypeName
      << " are "
      << msgList
      << "."
      << std::endl;
}
//--------------------------------------------------------------------------
void
SuppAlgBuilderLog::print_built_supp_alg_names(std::string suppAlgTypeName)
{
  const auto builtSuppAlgNames = built_supp_alg_names(suppAlgTypeName);

  std::string msgList = "";
  for (const auto& name : builtSuppAlgNames) {
    msgList += "`" + name + "', ";
  }
  msgList = msgList.substr(0, msgList.size()-2);

  NaluEnv::self().naluOutputP0()
      << "Built Supplemental Algortihms for "
      << suppAlgTypeName
      << " are "
      << msgList
      << "."
      << std::endl;
}
//--------------------------------------------------------------------------
std::set<std::string>
SuppAlgBuilderLog::valid_supp_alg_names(std::string suppAlgTypeName)
{
  auto it = validSuppAlgNames_.find(suppAlgTypeName);
  if (it == validSuppAlgNames_.end()) {
    return std::set<std::string>();
  }
  return it->second;
}
//--------------------------------------------------------------------------
std::set<std::string>
SuppAlgBuilderLog::built_supp_alg_names(std::string suppAlgTypeName)
{
  auto it = builtSuppAlgNames_.find(suppAlgTypeName);
  if (it == builtSuppAlgNames_.end()) {
    return std::set<std::string>();
  }
  return it->second;
}

} // namespace nalu
} // namespace Sierra
