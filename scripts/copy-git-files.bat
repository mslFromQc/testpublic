@echo off
setlocal enabledelayedexpansion
pushd %~dp0..

md compatibilities\Result

curl -o compatibilities/Result/Result.h https://srv-qc-rdgit.olympus-ossa.com/ondtlib/Utils/-/raw/master/src/Result/Result.h
curl -o compatibilities/Result/ResultDetails.h https://srv-qc-rdgit.olympus-ossa.com/ondtlib/Utils/-/raw/master/src/Result/ResultDetails.h
															      
popd
endlocal
@echo on
