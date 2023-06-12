#include <fstream>
#include <numbers>
#include <filesystem>

#include <hdf5/hdf5.h>
#include <rapidjson/document.h>

namespace rj = rapidjson;
namespace fs = std::filesystem;

bool WriteSetup(hid_t fileId_, const std::string& jsonSetup_);
std::vector<int16_t> CreateAscanAmplitudes(uint64_t dimSize0_, uint64_t dimSize1_, uint64_t dimSize2_);


int main()
{
  fs::path scenariosPathName = fs::current_path() / "../../scenarios";
  fs::path inputPathFilename = scenariosPathName / R"(general-mapping\oneline-scan-plate-paut-linear.json)";
  //fs::path inputPathFilename = scenariosPathName / R"(general-weld\oneline-scan-cod-paut-secorial.json)";
  
  fs::path outputPathFilename = fs::path(inputPathFilename).replace_extension("nde");
  hid_t fileId = H5Fcreate(outputPathFilename.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  std::ifstream ifs(inputPathFilename);
  std::string jsonSetup((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
  bool success = WriteSetup(fileId, jsonSetup);

  hid_t propId = H5Pcreate(H5P_LINK_CREATE);
  herr_t errorCode = H5Pset_create_intermediate_group(propId, 1);

  rj::Document jsonDoc;
  jsonDoc.Parse(jsonSetup.c_str());

  for (const auto& group : jsonDoc["groups"].GetArray())
  {
    if (group["dataset"].HasMember("ascan"))
    {
      // Create ascan amplitude dataset
      auto& amplitude = group["dataset"]["ascan"]["amplitude"];
      uint64_t dimSize0 = amplitude["dimensions"][0]["quantity"].GetUint64();

      uint64_t dimSize1{};
      std::string axisName = amplitude["dimensions"][1]["axis"].GetString();
      if (axisName == "VCoordinate") {
        dimSize1 = amplitude["dimensions"][1]["quantity"].GetUint64();
      }
      else if (axisName == "Beam") {
        dimSize1 = amplitude["dimensions"][1]["beams"].GetArray().Size();
      }

      uint64_t dimSize2 = amplitude["dimensions"][2]["quantity"].GetUint64();
      std::vector<uint64_t> ascanAmpDimSizes{ dimSize0, dimSize1, dimSize2 };

      hid_t amplitudeDspaceId = H5Screate_simple(static_cast<int>(ascanAmpDimSizes.size()), ascanAmpDimSizes.data(), nullptr);
      hid_t amplitudeDsetId = H5Dcreate(fileId, amplitude["path"].GetString(), H5T_NATIVE_SHORT, amplitudeDspaceId, propId, H5P_DEFAULT, H5P_DEFAULT);

      // Write amplitude dataset
      std::vector<int16_t> amplitudes = CreateAscanAmplitudes(dimSize0, dimSize1, dimSize2);

      errorCode = H5Dwrite(amplitudeDsetId, H5Dget_type(amplitudeDsetId), amplitudeDspaceId, amplitudeDspaceId, H5P_DEFAULT, amplitudes.data());
      errorCode = H5Sclose(amplitudeDspaceId);
      errorCode = H5Dclose(amplitudeDsetId);

      // Create ascan status dataset
      auto& status = group["dataset"]["ascan"]["status"];
      dimSize0 = status["dimensions"][0]["quantity"].GetUint64();

      axisName = status["dimensions"][1]["axis"].GetString();
      if (axisName == "VCoordinate") {
        dimSize1 = status["dimensions"][1]["quantity"].GetUint64();
      }
      else if (axisName == "Beam") {
        dimSize1 = status["dimensions"][1]["beams"].GetArray().Size();
      }

      std::vector<uint64_t> ascanStatusDimSizes{ dimSize0, dimSize1 };

      hid_t ascanStatusDspaceId = H5Screate_simple(static_cast<int>(ascanStatusDimSizes.size()), ascanStatusDimSizes.data(), nullptr);
      hid_t ascanStatusDsetId = H5Dcreate(fileId, status["path"].GetString(), H5T_NATIVE_B8, ascanStatusDspaceId, propId, H5P_DEFAULT, H5P_DEFAULT);

      // Write ascan status dataset
      std::vector<uint8_t> statuses(dimSize0 * dimSize1, 1);
      errorCode = H5Dwrite(ascanStatusDsetId, H5Dget_type(ascanStatusDsetId), ascanStatusDspaceId, ascanStatusDspaceId, H5P_DEFAULT, statuses.data());

      errorCode = H5Sclose(ascanStatusDspaceId);
      errorCode = H5Dclose(ascanStatusDsetId);
    }

    if (group["dataset"].HasMember("gateCscans"))
    {
      for (const auto& gateCscan : group["dataset"]["gateCscans"].GetArray())
      {
        // Create gate cscan dataset
        uint64_t dimSize0 = gateCscan["dimensions"][0]["quantity"].GetUint64();
        uint64_t dimSize1 = gateCscan["dimensions"][1]["quantity"].GetUint64();
        std::vector<uint64_t> gateCscanDimSizes{ dimSize0, dimSize1 };

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

        hid_t gateCscanDspaceId = H5Screate_simple(static_cast<int>(gateCscanDimSizes.size()), gateCscanDimSizes.data(), nullptr);
        hid_t gateCscanAmpDsetId = H5Dcreate(fileId, gateCscan["path"].GetString(), cscanDataTypeId, gateCscanDspaceId, propId, H5P_DEFAULT, H5P_DEFAULT);

        // Write gate cscan dataset
        std::vector<CompoundCscanData> compoundCscanData(dimSize0 * dimSize1);

        for (int32_t dim0{}; dim0 < dimSize0; dim0++) {
          for (int32_t dim1{}; dim1 < dimSize1; dim1++)
          {
            uint8_t status(1); // HasData
            int32_t peak = dim1 * 270;
            float peakTime(55.0e-6f + dim1 * 0.1e-6f);
            float crossingTime(40.0e-6f + dim1 * 0.1e-6f);

            auto cellIdx = dim1 + dim0 * dimSize1;
            compoundCscanData[cellIdx] = { crossingTime, peakTime, peak, status };
          }
        }

        errorCode = H5Dwrite(gateCscanAmpDsetId, H5Dget_type(gateCscanAmpDsetId), gateCscanDspaceId, gateCscanDspaceId, H5P_DEFAULT, compoundCscanData.data());

        errorCode = H5Tclose(cscanDataTypeId);
        errorCode = H5Sclose(gateCscanDspaceId);
        errorCode = H5Dclose(gateCscanAmpDsetId);
      }
    }
  }

  errorCode = H5Pclose(propId);
  errorCode = H5Fclose(fileId);

  return EXIT_SUCCESS;
}

std::vector<int16_t> CreateAscanAmplitudes(uint64_t dimSize0_, uint64_t dimSize1_, uint64_t dimSize2_)
{
  double repetition(2);
  double amplitude(_I16_MAX / 2);
  double offset(_I16_MAX / 2);
  std::vector<int16_t> amplitudes(dimSize0_ * dimSize1_ * dimSize2_);

  // Generate sinus waveforms retified
  for (size_t dim0{}; dim0 < dimSize0_; dim0++) {
    for (size_t dim1{}; dim1 < dimSize1_; dim1++)
    {
      auto index = dim1 * dimSize2_ + dim0 * dimSize1_ * dimSize2_;
      for (size_t dim2{}; dim2 < dimSize2_; dim2++)
      {
        double frequency = dim2 / (dimSize2_ / repetition - 1.0);
        double value = std::max(amplitude * std::sin(2 * std::numbers::pi * frequency), 0.0);
        amplitudes[index + dim2] = static_cast<int16_t>(std::lround(value));
      }
    }
  }

  return amplitudes;
}

bool WriteSetup(hid_t fileId_, const std::string& jsonSetup_)
{
  rj::Document schemaJsonDoc;
  fs::path rootPath = fs::current_path() / "../..";

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
