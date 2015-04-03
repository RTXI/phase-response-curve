/*
Copyright (C) 2011 Georgia Institute of Technology

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <scatterplot.h>
#include <runningstat.h>
#include <plotdialog.h>
//#include <RTXIprintfilter.h>
#include <default_gui_model.h>
#include <cstdlib>
#include <QtGui>
//using namespace std;

class PRC : public DefaultGUIModel {
	Q_OBJECT // macro needed if slots are implemented
	
	public:
		PRC(void);
		virtual ~PRC(void);
		void execute(void);
		void createGUI(DefaultGUIModel::variable_t *, int);
		void customizeGUI(void);
		
		enum mode_t {
			PRC_BASELINE, // measuring intrinsic P0
			PRC_PRE, // beginning of protocol, count 10 periods
			PRC_POST, // after perturbed period
			PRC_DONE,
		};
		
		enum extra_flags_t {
			STARTBASELINE,
			MEASUREPRC,
		};
	
	public slots:
	
	signals:
		void setPlotRange(double xmin, double xmax, double ymin, double ymax);
		void saveOK(bool);
		void setPlotMode(bool);
		void drawData(double* x, double* y, int size);
		void drawLine(double* x, double* y, int size);
	
	protected:
		virtual void update(DefaultGUIModel::update_flags_t);
		void update(PRC::extra_flags_t);
	
	private:
		// inputs, states, calculated values
		double systime;
		double Vm;
		double dt;
		// parameters
		double mindelay;
		double maxdelay;
		double newmindelay;
		double newmaxdelay;
		double stepsize; // step minphase and maxphase
		double gmax; // maximum conductance for alpha shaped perturbation
		double tau; // time constant
		double esyn; // synaptic reversal potential
		double repeat;
		// options
		bool randomize;
		bool plotPRC;
		bool saveData;
		mode_t mode;
		//bookkeeping
		double spktime;
		double prevspktime;
		int nstep;
		double ISI; // = spktime - prevspktime
		double P0; // intrinsic period
		double P1; // period of perturbed cycle
		double P2; // period of cycle following perturbed cycle
		bool random; // is shuffled?
		double periodnum;
		int Gsyncount;
		int stepcount; // which step within cycle
		int cyclecount; // number of repeated cycles
		long long count; // keep track of time
		int spikecount;
		int nidx; // arbitrary idx
		
		RunningStat runningPeriod;
		
		double* arrGsyn; // holds the stimulus conductance values
		double* arrISI; // holds five periods for averaging
		std::vector<double> reftr; // hold values for reference tr-ts plots
		std::vector<double> refts;
		
		double* delays; // holds delays for each protocol
		// these hold the data to be plotted until it is manually cleared
		std::vector<double> arrtr; // milliseconds
		std::vector<double> arrts; // holds delays in order for each cycle in seconds
		std::vector<double> arrP0; // milliseconds
		std::vector<double> arrphi; // holds phases, unitless
		std::vector<double> arrPRC1; // store 1st order PRC, unitless
		std::vector<double> arrPRC2; // store 2nd order PRC, unitless
		
		// Spike Detector parameters for cell
		double thresh;
		double min_int;
		int spikestate;
		double last_spike;
		
		// PRC functions
		void initDelayArray(); // initialize stimulus delays
		void initParameters(); // initialize parameters
		void initStimulus(); // initialize stimulus conductance waveform based on parameters
		void bookkeep(); // reset bookkeeping variables
		void countspikes(); // increment spikes, cycles, compute ISI
		void resetcycle(); // bookkeeping
		void shuffleonce(); // shuffle stimulus delays once per cycle
		double shiftstimulus(double spiketime); // determine correct output amplitude
		void makeplot();
		//    void drawRefLines();
		//    void LoadRefData(QString);
		
		// QT components
		ScatterPlot *splot;
		QLabel *moduleStatus;
		QButtonGroup *plotGroup;
		
		// Functions and parameters for saving data to file without using data recorder
		bool OpenFile(QString);
		QFile dataFile;
		QTextStream stream;
	
	private slots:
		// custom
		void togglerandom(bool);
		void togglePlot(bool);
		void savePRC();
		void startBaseline();
		void measurePRC();
		void print();
		void exportSVG();
		void clearData();
};
