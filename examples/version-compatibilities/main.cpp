#include <iostream>
#include <filesystem>

#include <hdf5/hdf5.h>
#include <rapidjson/schema.h>
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>

#include "compatibilities/FileFormat.h"
#include "compatibilities/FileFormatChanges.h"

namespace rj = rapidjson;
namespace fs = std::filesystem;


int main()
{
  std::cout << "Enter the path filename: " << std::endl;

  std::string inputPathFilename;
  std::getline(std::cin, inputPathFilename);
  inputPathFilename.erase(remove(inputPathFilename.begin(), inputPathFilename.end(), '\"'), inputPathFilename.end());
  if (!fs::exists(inputPathFilename)) {
    return EXIT_FAILURE;
  }

  // Open file and JSON dataset
  hid_t fileId = H5Fopen(inputPathFilename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t propId = H5Pcreate(H5P_DATASET_ACCESS);
  hid_t dsetId = H5Dopen(fileId, "/Domain/Setup", propId);

  // Read JSON metadata
  hid_t dataType = H5Dget_type(dsetId);
  std::vector<char> json(H5Tget_size(dataType));
  herr_t errorCode = H5Dread(dsetId, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, json.data());
  std::string jsonSetupFile(json.begin(), json.end());
  errorCode = H5Dclose(dsetId);
  errorCode = H5Fclose(fileId);


  std::string errorMessage;
  fs::path rootPath = fs::current_path() / "../..";
  fs::path schemasPath = rootPath / "schemas";

  NDE::TResult<std::string> jsonSetup = NDE::Upgrade(schemasPath, jsonSetupFile);
  if (jsonSetup.IsSuccess())
  {
    fs::path jsonPathFilename(inputPathFilename);
    std::ofstream out(jsonPathFilename.parent_path() / "setup-updated.json");
    out << jsonSetup.Value();
    out.close();
  }
  else {
    std::cout << errorMessage;
  }

  return EXIT_SUCCESS;
}
