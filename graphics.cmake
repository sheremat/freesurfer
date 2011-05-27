FIND_PACKAGE(ITK)
IF(ITK_FOUND)
   INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

FIND_PACKAGE(VTK)
IF(VTK_FOUND) 
   INCLUDE(${VTK_USE_FILE})
ENDIF(VTK_FOUND)

FIND_PACKAGE(KWWidgets)
IF(KWWidgets_FOUND)
  INCLUDE(${KWWidgets_USE_FILE})
ENDIF(KWWidgets_FOUND)

FIND_PACKAGE (TCL)
FIND_PACKAGE (OPENGL)

SET (VTK_LIBS
  VTKCommon
)

SET (ITK_LIBS
  ITKCommon
  ITKIO
)

SET (GRAPH_INCLUDES
  ${VTK_INCLUDE_DIRS}
  ${KWWidgets_INCLUDE_DIRS}
  ${ITK_INCLUDE_DIRS}
  ${TCL_INCLUDE_PATH}
  ${OPENGL_INCLUDE_DIR}
)
