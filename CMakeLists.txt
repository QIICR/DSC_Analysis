cmake_minimum_required(VERSION 2.8.7)

if (NOT Slicer_SOURCE_DIR)
  set(EXTENSION_NAME DSCMRIAnalysis)
  set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Modules/DSC_MRI_Analysis")
  set(EXTENSION_CATEGORY "Quantification")
  set(EXTENSION_CONTRIBUTORS "Xiao Da (MGH)")
  set(EXTENSION_DESCRIPTION "Analysis of Dynamic Susceptibility Contrast (DSC) MRI.")
  set(EXTENSION_ICONURL "http://slicer.org/slicerWiki/images/9/92/DSC_logo_Resized.png")
  set(EXTENSION_STATUS "Beta")
  set(EXTENSION_SCREENSHOTURLS "http://slicer.org/slicerWiki/images/e/e8/DSC_GUI_V0.png")
  set(EXTENSION_ACKNOWLEDGEMENTS "This work is part of the National Alliance for Medical Image Computing (NA-MIC), funded by the National Institutes of Health through the NIH Roadmap for Medical Research, and by National Cancer Institute as part of the Quantitative Imaging Network initiative (U01CA154601) and QIICR (U24CA180918).")
endif()

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()

#-----------------------------------------------------------------------------
# PkModeling has a support library that needs to be installed as part
# of the packaging
#
# Disabling.  Changed the PkSolver library to be build statiically so it does not need to be installed
# option(PkModeling_BUILD_SLICER_EXTENSION "Build as a Slicer extension" OFF)
#
# set(PkModeling_INSTALL_SUPPORT_RUNTIME_DESTINATION bin)
# set(PkModeling_INSTALL_SUPPORT_LIBRARY_DESTINATION lib)
# set(PkModeling_INSTALL_SUPPORT_ARCHIVE_DESTINATION lib)
# if(PkModeling_BUILD_SLICER_EXTENSION)
#  set(PkModeling_INSTALL_SUPPORT_RUNTIME_DESTINATION ${Slicer_INSTALL_BIN_DIR})
#  set(PkModeling_INSTALL_SUPPORT_LIBRARY_DESTINATION ${Slicer_INSTALL_LIB_DIR})
#  set(PkModeling_INSTALL_SUPPORT_ARCHIVE_DESTINATION ${Slicer_INSTALL_LIB_DIR})
# endif()

#-----------------------------------------------------------------------------
add_subdirectory(PkSolver)
add_subdirectory(CLI)

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()

