
#ifndef VANDLE_MERGER_YAMLPARAMETER_H
#define VANDLE_MERGER_YAMLPARAMETER_H

#include "yaml-cpp/yaml.h"

class YamlParameter : public YAML::Node{
private:
  YamlParameter(std::string file_name) : YAML::Node( YAML::LoadFile(file_name) ) {}
  ~YamlParameter() = default; 

  static YamlParameter* instance;

public:
  YamlParameter() = delete;
  YamlParameter(const YamlParameter&) = delete;
  YamlParameter& operator = (const YamlParameter&) = delete;
  YamlParameter(YamlParameter&&) = delete;
  YamlParameter& operator = (YamlParameter&&) = delete;

  static YamlParameter& GetInstance() {
    return *instance;
  }

  static void Create(std::string file_name) {
    if( !instance ) {
      instance = new YamlParameter(file_name);
    }
  }

  static void Destroy() {
    if( instance ) {
      delete instance;
      instance = NULL;
    }
  }
};

#endif /** VANDLE_MERGER_YAMLPARAMETER_H **/
