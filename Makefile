PLUGIN_NAME = phase_response_curve

RTXI_INCLUDES=/usr/local/lib/rtxi_includes

HEADERS = phase-response-curve.h\
          $(RTXI_INCLUDES)/runningstat.h\
          $(RTXI_INCLUDES)/basicplot.h\
          $(RTXI_INCLUDES)/scatterplot.h\
          $(RTXI_INCLUDES)/incrementalplot.h\
          $(RTXI_INCLUDES)/scrollbar.h\
          $(RTXI_INCLUDES)/scrollzoomer.h\

SOURCES = phase-response-curve.cpp \
          moc_phase-response-curve.cpp\
          $(RTXI_INCLUDES)/runningstat.cpp\
          $(RTXI_INCLUDES)/basicplot.cpp\
          $(RTXI_INCLUDES)/scatterplot.cpp\
          $(RTXI_INCLUDES)/incrementalplot.cpp\
          $(RTXI_INCLUDES)/scrollbar.cpp\
          $(RTXI_INCLUDES)/scrollzoomer.cpp\
          $(RTXI_INCLUDES)/moc_scatterplot.cpp\
          $(RTXI_INCLUDES)/moc_incrementalplot.cpp\
          $(RTXI_INCLUDES)/moc_scrollbar.cpp\
          $(RTXI_INCLUDES)/moc_scrollzoomer.cpp\
          $(RTXI_INCLUDES)/moc_basicplot.cpp

LIBS = -lqwt

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile
