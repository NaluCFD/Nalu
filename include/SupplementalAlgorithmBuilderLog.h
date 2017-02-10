/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#ifndef SupplementalAlgorithmBuilderLog_h
#define SupplementalAlgorithmBuilderLog_h

#include <map>
#include <set>
#include <string>
#include <vector>

namespace sierra{
namespace nalu{
  class SuppAlgBuilderLog
  {
  public:
    static SuppAlgBuilderLog& self();
    SuppAlgBuilderLog(const SuppAlgBuilderLog&) = delete;
    void operator=(const SuppAlgBuilderLog&) = delete;

    void add_valid_name(std::string suppAlgTypeName, std::string name);
    void add_built_name(std::string suppAlgTypeName, std::string name);

    bool print_invalid_supp_alg_names(
      std::string suppAlgTypeName,
      const std::map<std::string, std::vector<std::string>>& inputFileNames );

    void print_valid_supp_alg_names(std::string suppAlgTypeName);
    void print_built_supp_alg_names(std::string suppAlgTypeName);

    std::set<std::string> valid_supp_alg_names(std::string suppAlgTypeName);
    std::set<std::string> built_supp_alg_names(std::string suppAlgTypeName);
  private:
    SuppAlgBuilderLog() = default;

    std::map<std::string, std::set<std::string>> validSuppAlgNames_;
    std::map<std::string, std::set<std::string>> builtSuppAlgNames_;
  };

} // namespace nalu
} // namespace Sierra

#endif
