PLUGIN_NAME = PRC

HEADERS = PRC.h\
          include/runningstat.h\
          include/scrollzoomer.h\
          include/scatterplot.h\
          include/basicplot.h\
          include/incrementalplot.h\
          include/scrollbar.h\

SOURCES = PRC.cpp \
          moc_PRC.cpp\
          include/runningstat.cpp\
          include/basicplot.cpp\
          include/scatterplot.cpp\
          include/incrementalplot.cpp\
          include/scrollbar.cpp\
          include/scrollzoomer.cpp\
          include/moc_runningstat.cpp\
          include/moc_scatterplot.cpp\
          include/moc_incrementalplot.cpp\
          include/moc_scrollbar.cpp\
          include/moc_scrollzoomer.cpp\
          include/moc_basicplot.cpp

LIBS = -lqwt

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile
