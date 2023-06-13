#include <fstream>
#include <numbers>
#include <filesystem>

#include <hdf5/hdf5.h>
#include <rapidjson/schema.h>
#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>

namespace rj = rapidjson;
namespace fs = std::filesystem;

bool WriteSetup(hid_t fileId_, const std::string& jsonSetup_);


int main()
{
  fs::path scenariosPathName = fs::current_path() / "../../scenarios";
  fs::path jsonPathFilename = scenariosPathName / "general-mapping/gate-cscan-only.json";
  std::ifstream ifs(jsonPathFilename);
  std::string jsonSetup((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));

  fs::path outputPathFilename = jsonPathFilename.replace_extension("nde");
  hid_t fileId = H5Fcreate(outputPathFilename.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  hid_t propId = H5Pcreate(H5P_LINK_CREATE);
  herr_t errorCode = H5Pset_create_intermediate_group(propId, 1);

  struct CompoundCscanData
  {
    float CrossingTime{};
    float PeakTime{};
    int Peak{};
    uint8_t Status{};
  };

  hid_t cscanDataTypeId = H5Tcreate(H5T_COMPOUND, sizeof(CompoundCscanData));
  H5Tinsert(cscanDataTypeId, "crossingTime", HOFFSET(CompoundCscanData, CrossingTime), H5T_NATIVE_FLOAT);
  H5Tinsert(cscanDataTypeId, "peakTime", HOFFSET(CompoundCscanData, PeakTime), H5T_NATIVE_FLOAT);
  H5Tinsert(cscanDataTypeId, "peak", HOFFSET(CompoundCscanData, Peak), H5T_NATIVE_INT);
  H5Tinsert(cscanDataTypeId, "status", HOFFSET(CompoundCscanData, Status), H5T_NATIVE_B8);

  rj::Document jsonDoc;
  jsonDoc.Parse(jsonSetup.c_str());
  bool success = WriteSetup(fileId, jsonSetup);

  for (const auto& group : jsonDoc["groups"].GetArray())
  {
    for (const auto& gateCscan : group["dataset"]["gateCscans"].GetArray())
    {
      // Create dataset
      uint64_t dimSize0 = gateCscan["dimensions"][0]["quantity"].GetUint64();
      uint64_t dimSize1 = gateCscan["dimensions"][1]["quantity"].GetUint64();
      std::vector<uint64_t> cscanDimSizes{ dimSize0, dimSize1 };

      hid_t cscanDspaceId = H5Screate_simple(static_cast<int>(cscanDimSizes.size()), cscanDimSizes.data(), nullptr);
      hid_t cscanDsetId = H5Dcreate(fileId, gateCscan["path"].GetString(), cscanDataTypeId, cscanDspaceId, propId, H5P_DEFAULT, H5P_DEFAULT);

      // Write dataset
      std::vector<CompoundCscanData> compoundCscanData(dimSize0 * dimSize1);

      for (int32_t dim0{}; dim0 < dimSize0; dim0++) {
        for (int32_t dim1{}; dim1 < dimSize1; dim1++)
        {
          uint8_t status(1); // HasData
          int32_t peak = _I16_MIN + dim1 * 1024;
          float peakTime(55.0e-6f + dim1 * 0.1e-6f);
          float crossingTime(40.0e-6f + dim1 * 0.1e-6f);

          auto cellIdx = dim1 + dim0 * dimSize1;
          compoundCscanData[cellIdx] = { crossingTime, peakTime, peak, status };
        }
      }

      errorCode = H5Dwrite(cscanDsetId, H5Dget_type(cscanDsetId), cscanDspaceId, cscanDspaceId, H5P_DEFAULT, compoundCscanData.data());
      errorCode = H5Pclose(propId);
      errorCode = H5Tclose(cscanDataTypeId);
      errorCode = H5Dclose(cscanDsetId);
    }
  }

  errorCode = H5Fclose(fileId);

  return EXIT_SUCCESS;
}

bool WriteSetup(hid_t fileId_, const std::string& jsonSetup_)
{
  rj::Document setupJsonDoc;
  if (!setupJsonDoc.Parse(jsonSetup_.c_str()).HasParseError())
  {
    // Write JSON dataset
    hid_t dtype = H5Tcopy(H5T_C_S1);
    herr_t errorCode = H5Tset_size(dtype, jsonSetup_.length());
    errorCode = H5Tset_cset(dtype, H5T_CSET_UTF8);
    errorCode = H5Tset_strpad(dtype, H5T_STR_NULLTERM);

    hid_t dspaceId = H5Screate(H5S_SCALAR);
    hid_t propId = H5Pcreate(H5P_LINK_CREATE);
    errorCode = H5Pset_create_intermediate_group(propId, 1);
    hid_t dsetId = H5Dcreate(fileId_, "/Domain/Setup", dtype, dspaceId, propId, H5P_DEFAULT, H5P_DEFAULT);
    H5Pclose(propId);

    errorCode = H5Dwrite(dsetId, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, jsonSetup_.data());
    H5Dclose(dsetId);

    // Write version attribute
    std::string version = setupJsonDoc["version"].GetString();
    errorCode = H5Tset_size(dtype, version.length());
    hid_t attributeId = H5Acreate(fileId_, "Format version", dtype, dspaceId, H5P_DEFAULT, H5P_DEFAULT);
    errorCode = H5Awrite(attributeId, dtype, version.c_str());
    H5Aclose(attributeId);
    H5Sclose(dspaceId);

    return true;
  }

  return false;
}
