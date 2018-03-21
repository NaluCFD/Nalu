/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#ifndef KernelBuilderLog_h
#define KernelBuilderLog_h

#include <map>
#include <set>
#include <string>
#include <vector>

namespace sierra{
namespace nalu{
  class KernelBuilderLog
  {
  public:
    static KernelBuilderLog& self();
    KernelBuilderLog(const KernelBuilderLog&) = delete;
    void operator=(const KernelBuilderLog&) = delete;

    void add_valid_name(std::string kernelTypeName, std::string name);
    void add_built_name(std::string kernelTypeName, std::string name);

    bool print_invalid_kernel_names(
      std::string kernelTypeName,
      const std::map<std::string, std::vector<std::string>>& inputFileNames );

    void print_valid_kernel_names(std::string kernelTypeName);
    void print_built_kernel_names(std::string kernelTypeName);

    std::set<std::string> valid_kernel_names(std::string kernelTypeName);
    std::set<std::string> built_kernel_names(std::string kernelTypeName);
  private:
    KernelBuilderLog() = default;

    std::map<std::string, std::set<std::string>> validKernelNames_;
    std::map<std::string, std::set<std::string>> builtKernelNames_;
  };

} // namespace nalu
} // namespace Sierra

#endif
