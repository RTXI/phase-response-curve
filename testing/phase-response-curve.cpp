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

#include "phase-response-curve.h"
#include <math.h>
#include <algorithm>
#include <time.h>

#include <QtGui>
#include <sys/stat.h>

#if QT_VERSION >= 0x040300
#ifdef QT_SVG_LIB
#include <QSvgGenerator>
#endif
#endif
#if QT_VERSION >= 0x040000
#include <QPrintDialog>
#include <QFileInfo>
#else
#include <qwt_painter.h>
#endif
//#include <qwt_array.h>

extern "C" Plugin::Object *createRTXIPlugin(void) {
	return new PRC();
}

// inputs, outputs, parameters for default_gui_model
static DefaultGUIModel::variable_t vars[] = {
	{ "ISI (ms)", "ISI (ms)", DefaultGUIModel::STATE, },
	{ "Period Number", "Period Number", DefaultGUIModel::STATE, },
	{ "Vm", "Membrane voltage of real cell (V)", DefaultGUIModel::INPUT, },
	{ "Command", "Command perturbation current", DefaultGUIModel::OUTPUT, },
	{ "Cell spike state", "Real cell spike state", DefaultGUIModel::OUTPUT, },
	{ "Min Delay (ms)", "Minimum stimulus delay (ms)",
		DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Max Delay (ms)", "Maximum stimlus delay (ms)",
		DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Step Size (ms)", "Step size between minimum and maximum delay (ms)",
		DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Gmax (nS)", "Maximum synaptic conductance for stimulus",
		DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Time Constant tau (ms)", "Time constant for alpha-shaped conductance",
		DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Esyn (mV)", "Reversal potential for stimulus",
		DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Repeat", "Number of times to run cycle", DefaultGUIModel::PARAMETER
		| DefaultGUIModel::DOUBLE, },
	{ "Threshold (mV)", "Threshold (mV) at which to detect a spike",
		DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Min Interval (s)", "Minimum interval(s) that must pass between spikes",
		DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE, },
	{ "Time (s)", "Time", DefaultGUIModel::STATE, },
};

// some necessary variable DO NOT EDIT
static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

// Default constructor
PRC::PRC(void) : DefaultGUIModel("Phase Response Curve", ::vars, ::num_vars) {
	setWhatsThis(
	"<p><b>Phase Response Curve:</b></p><p> This module applies an alpha-shaped conductance to "
	"the cell at a fixed delay. For the PRC protocol, this occurs after specified number of unperturbed "
	"interspike intervals (ISI), during which an intrinsic period P0 is computed by averaging the ISIs. "
	"The period after the perturbed period is designated as P1 and the period following P1 is designated "
	"as P2. The first and second order PRCs are computed as PRC1=(P1-P0)/P0 and PRC2=(P2-P0)/P0. From this "
	"data, the module plots the PRC curves in phase units or as ts-tr curves (stimulus time - response time)."
	"<br><br>PRC Functions:<br>Measure intrinsic P0: computes running average of the measured ISIs to help you choose "
	"the fixed delays at which to apply the perturbation.<br>"
	"Measure PRC: Runs the PRC protocol defined by the parameters.<br><br>Input(0): Vm in (V) of the cell "
	"whose PRC you want to measure.<br><br>Output(0): PRC command current</p>");

	// initialize
	initParameters();
	initDelayArray();
	initStimulus();
	DefaultGUIModel::createGUI(vars, num_vars);
	customizeGUI();
	update( INIT );
	refresh();
	QTimer::singleShot(0, this, SLOT(resizeMe()));
	printf("\nPRC Module loaded...\n");
}

PRC::~PRC(void) {}

void PRC::execute(void) {
	Vm = input(0); // input in V, real cell's Vm
	systime = count * dt; // current time
	
	// handle spike detection
	switch (spikestate) {
		case 0:
			if (Vm > thresh){
				spikestate = 1;
				last_spike = systime;
			}
			break;
		
		case 1:
			spikestate = 2;
			break;
	
		case 2:
			if (Vm > thresh && (systime - last_spike) > 100) spikestate = 4;
			else if (Vm < thresh) spikestate = 3;
			break;

		case 3:
			spikestate = -1;
			break;

		case 4:
			if (Vm < thresh) spikestate = -1; 
			break;

		case -1:
			if (systime - last_spike > min_int) spikestate = 0;
			break;

		default:
			break;
	}
	output(1) = spikestate;
	
	// compute current to output
	switch (mode) {
		case PRC_BASELINE: // get intrinsic period
			if (spikestate == 1) countspikes();
			if (spikecount > 2) runningPeriod.push(ISI); // ignore first period}
			output(0) = 0;
			break; // end case mode=PRC_BASELINE

		case PRC_PRE: // count 10 periods, calculate first P0, switch mode
			if (spikestate == 1) {
				countspikes();
				if (spikecount >= 7 && spikestate == 1) runningPeriod.push(ISI); // save last 5 periods (#6-10) to get P0
				if (spikecount == 11) { // period #10
					P0 = runningPeriod.mean();
					arrP0[stepcount] = P0;
					arrts[stepcount] = delays[stepcount];
					arrphi[stepcount] = arrts[stepcount] / P0;
					printf("-P0: %f, ts: %f, phi: %f\n", P0, arrts[stepcount],
					arrphi[stepcount]);
					runningPeriod.clear();
					output(0) = shiftstimulus(spktime);
					mode = PRC_POST;
				}
			}
			output(0) = 0;
			break; // end case mode=PRC_PRE

		case PRC_POST:
			if (spikestate == 1) { // if spike just occurred, count spikes, compute PRCs
				countspikes();
				if (cyclecount < repeat) { // as long as there are cycles remaining
					shuffleonce();
					if (spikecount >= 7) // save only last 5 periods (#6-10) for averaging to get P0
						runningPeriod.push(ISI);
					//printf("Period Number: %i, ISI: %f, ISIavg: %f\n", spikecount - 1, ISI, runningPeriod.mean());
					switch (spikecount) { // calculate PRC and do bookkeeping
						case 11: // spike #11, currently period #11, perturbed cycle, get P0
							Gsyncount = 0;
							P0 = runningPeriod.mean();
							arrP0[stepcount] = P0;
							arrts[stepcount] = delays[stepcount];
							arrphi[stepcount] = arrts[stepcount] / P0;
							printf("-P0: %f, ts: %f, phi: %f\n", P0, arrts[stepcount],
							arrphi[stepcount]);
							break;

						case 12: // spike #12, currently period #12, cycle after perturbed cycle, get P1
							P1 = ISI;
							nidx = static_cast<int> (cyclecount * nstep + stepcount);
							arrPRC1[stepcount] = (P1 - P0) / P0;
							arrtr[stepcount] = ISI - arrts[stepcount];
							printf(" P1: %f, tr: %f, F1: %f\n", P1, arrtr[stepcount],
							arrPRC1[stepcount]);
							spikecount = 1; // reset spike count for next stimulus delay
							runningPeriod.clear();
							break;

						case 2: // spike #13, currently period #13, 2 cycles after perturb cycle, get P2
							P2 = ISI;
							nidx = static_cast<int> (cyclecount * nstep + stepcount);
							arrPRC2[stepcount] = (P2 - P0) / P0;
							printf(" P2: %f, F2: %f\n", P2, arrPRC2[stepcount]);
							stepcount++;
							moduleStatus->setText("Measuring PRC... Delay = " + QString::number(delays[stepcount] / 1000) + "s");
							if (stepcount == nstep) resetcycle();
							break;
					
						default:
							break;
						} // end switch(spikecount)
				}
				else { // all cycles are done
					mode = PRC_DONE; // last PRC has been calculated
				}
			} // end if (spikestate == 1)

			switch (spikecount) { // between spikes, calculate correct output
				case 11: // spike #11, currently period #11, satisfied by first spike in PRC_POST, perturbed cycle, get P0
					output(0) = shiftstimulus(spktime);
					break;
			
				case 12: // spike #12, currently period #12, cycle after perturbed cycle, get P1, stimulus may extend into this period
					output(0) = shiftstimulus(prevspktime);
					break;
			
				case 2: // spike #13, currently period #13, 2 cycles after perturbed cycle, get P2
			
				default:
					output(0) = 0;
					break;
			} // end switch(spikecount)
			break; // end case mode=PRC_POST
		
		case PRC_DONE: // spike #13, currently period #12, cycle after perturbed cycle, get P1
			setActive(false);
			update( PAUSE );
			if (plotPRC == true) makeplot();
			break; // end case mode=PRC_DONE

		default:
			break;
	} // end switch(mode)
	
	count++; // increment time
	return;
} // end PRC::execute

void PRC::update(DefaultGUIModel::update_flags_t flag) {
	switch (flag) {
		case INIT:
			setState("ISI (ms)", ISI);
			setState("Period Number", periodnum);
			setParameter("Min Delay (ms)", QString::number(mindelay * 1e3)); // initialized in s, display in ms
			setParameter("Max Delay (ms)", QString::number(maxdelay * 1e3)); // initialized in s, display in ms
			setParameter("Step Size (ms)", QString::number(stepsize * 1e3)); // initialized in s, display in ms
			setParameter("Gmax (nS)", QString::number(gmax * 1e9)); // initialized in S, display in nS
			setParameter("Time Constant tau (ms)", QString::number(tau * 1e3)); // initialized in s, display in ms
			setParameter("Esyn (mV)", QString::number(esyn * 1e3)); // initialized in V, display in mV
			setParameter("Repeat", QString::number(repeat)); // initially 1
			setState("Time (s)", systime);
			
			// spike detector stuff
			setParameter("Threshold (mV)", QString::number(thresh * 1000.0)); // stored as V, display in mV
			setParameter("Min Interval (s)", QString::number(min_int));
			emit PRC::setPlotRange(0, 1, -1.1, 1.1);
			break;

		case MODIFY:
			if (getParameter("Time Constant tau (ms)").toDouble() <= 0) {
				QMessageBox::critical(this, "PRC: Modifying tau",
				"The time constant must have a positive value.\n");
				setParameter("Time Constant tau (ms)", QString::number(tau * 1e3));
			} else {
				tau = getParameter("Time Constant tau (ms)").toDouble() * 1e-3; // set by user in ms, change to s
			}
			mindelay = getParameter("Min Delay (ms)").toDouble() * 1e-3; // set by user in ms, change to s
			maxdelay = getParameter("Max Delay (ms)").toDouble() * 1e-3; // set by user in ms, change to s
			stepsize = getParameter("Step Size (ms)").toDouble() * 1e-3; // set by user in ms, change to s
			gmax = getParameter("Gmax (nS)").toDouble() * 1e-9; // set by user in nS, change to S
			esyn = getParameter("Esyn (mV)").toDouble() * 1e-3; // set by user in mV, change to V
			repeat = getParameter("Repeat").toDouble();
			
			// spike detector stuff
			thresh = getParameter("Threshold (mV)").toDouble() / 1000; // specified in mV, convert to V
			min_int = getParameter("Min Interval (s)").toDouble();
			bookkeep();
			nstep = int((maxdelay + 1e-12 - mindelay) / stepsize) + 1; // recalculate
			initDelayArray();
			initStimulus();
			break;

		case PAUSE:
			output(0) = 0; // stop command in case pause occurs in the middle of command
			switch (mode) {
				case PRC_BASELINE:
					moduleStatus->setText("Done getting intrinsic period.");
					printf("%f ms over %i periods\n", runningPeriod.mean(), spikecount - 1);
					break;
			
				case PRC_DONE:
					moduleStatus->setText("Done measuring PRC.");
					break;

				default:
					moduleStatus->setText("Ready.");
					break;
			}
			break;

		case UNPAUSE:
			bookkeep();
			initDelayArray();
			switch (mode) {
				case PRC_BASELINE:
					moduleStatus->setText("Getting intrinsic period...");
					break;
		
				case PRC_PRE:
					printf("Starting PRC protocol: mindelay: %f maxdelay: %f stepsize %f gmax %f tau %f esyn %f\n", mindelay * 1e3, maxdelay * 1e3, stepsize * 1e3, gmax * 1e-12, tau * 1e3, esyn * 1e3);
					break;

				case PRC_POST:
					moduleStatus->setText("Measuring PRC... Delay = " + 
					                      QString::number(delays[stepcount] / 1000) + "s");
					break;
			
				default:
					moduleStatus->setText("something is running, but what?");
					break;
			}
			break;

		case PERIOD:
			dt = RT::System::getInstance()->getPeriod() * 1e-9; // period in s
			initStimulus();
			break;

		default:
			break;
	}
}

void PRC::update(PRC::extra_flags_t flag) {
	switch (flag) {
		case STARTBASELINE:
			mode = PRC_BASELINE;
			bookkeep();
			printf("Getting intrinsic period: ");
			pause(false);
			moduleStatus->setText("Getting intrinsic period...");
			break;

		case MEASUREPRC:
			mode = PRC_PRE;
			bookkeep();
			printf("Measuring PRC: ");
			pause(false);
			moduleStatus->setText("Measuring PRC...");
			break;
	}
}

void PRC::customizeGUI(void) {
	QGridLayout *customLayout = DefaultGUIModel::getLayout();

	QGroupBox *plotBox = new QGroupBox("PRC Plot");
	QHBoxLayout *plotBoxLayout = new QHBoxLayout;
	plotBox->setLayout(plotBoxLayout);
	QPushButton *clearButton = new QPushButton("&Clear Data");
	QPushButton *savePlotButton = new QPushButton("Save Screenshot");
	QPushButton *printButton = new QPushButton("Print");
	QPushButton *saveDataButton = new QPushButton("&Save PRC Data");
	plotBoxLayout->addWidget(clearButton);
	plotBoxLayout->addWidget(savePlotButton);
	plotBoxLayout->addWidget(printButton);
	plotBoxLayout->addWidget(saveDataButton);

	splot = new ScatterPlot(this);
	QObject::connect(clearButton, SIGNAL(clicked()), splot, SLOT(clear()));
	QObject::connect(clearButton, SIGNAL(clicked()), this, SLOT(clearData()));
	QObject::connect(savePlotButton, SIGNAL(clicked()), this, SLOT(exportSVG()));
	QObject::connect(printButton, SIGNAL(clicked()), this, SLOT(print()));
	QObject::connect(saveDataButton, SIGNAL(clicked()), this, SLOT(savePRC()));

	clearButton->setToolTip("Clear");
	savePlotButton->setToolTip("Save screenshot");
	saveDataButton->setToolTip("Save PRC data to an ASCII file");
	printButton->setToolTip("Print plot");

	QObject::connect(this, SIGNAL(drawData(double*,double*,int)), splot, SLOT(appendData(double*,double*,int)));

	QGroupBox *baselineBox = new QGroupBox("PRC Functions");
	QHBoxLayout *baselineBoxLayout = new QHBoxLayout;
	baselineBox->setLayout(baselineBoxLayout);
	QButtonGroup *baselineButtons = new QButtonGroup;
	baselineButtons->setExclusive(true);
	QPushButton *startBttn = new QPushButton("Measure intrinsic P0");
	QPushButton *getPRCBttn = new QPushButton("Measure PRC");
	baselineBoxLayout->addWidget(startBttn);
	baselineBoxLayout->addWidget(getPRCBttn);
	baselineButtons->addButton(startBttn);
	baselineButtons->addButton(getPRCBttn);
	startBttn->setToolTip("Start running average of ISIs");
	getPRCBttn->setToolTip("Run PRC protocol");
	QObject::connect(startBttn, SIGNAL(clicked()), this, SLOT(startBaseline()));
	QObject::connect(getPRCBttn, SIGNAL(clicked()), this, SLOT(measurePRC()));
	
	// add custom components to layout below default_gui_model components
	QVBoxLayout *optionRowLayout = new QVBoxLayout;
	QCheckBox *randomCheckBox = new QCheckBox("Randomize");
	optionRowLayout->addWidget(randomCheckBox);
	QObject::connect(randomCheckBox, SIGNAL(toggled(bool)), this, SLOT(togglerandom(bool)));
	randomCheckBox->setToolTip("Randomize order of presentation");
	
	QGroupBox *plotTypeBox = new QGroupBox("Plot");
	QHBoxLayout *plotTypeBoxLayout = new QHBoxLayout;
	plotTypeBox->setLayout(plotTypeBoxLayout);
	plotGroup = new QButtonGroup;
	plotGroup->setExclusive(true);
	QRadioButton *plotNoneRadio = new QRadioButton("None");
	QRadioButton *plotPRCRadio = new QRadioButton("PRC");
	QRadioButton *plotTRTSRadio = new QRadioButton("ts-tr");
	plotTypeBoxLayout->addWidget(plotNoneRadio);
	plotTypeBoxLayout->addWidget(plotPRCRadio);
	plotTypeBoxLayout->addWidget(plotTRTSRadio);
	plotGroup->addButton(plotNoneRadio);
	plotGroup->addButton(plotPRCRadio);
	plotGroup->addButton(plotTRTSRadio);

	QObject::connect(plotPRCRadio, SIGNAL(toggled(bool)), splot, SLOT(setVisible(bool)));
	QObject::connect(plotPRCRadio, SIGNAL(toggled(bool)), plotBox, SLOT(setVisible(bool)));
	QObject::connect(plotPRCRadio, SIGNAL(toggled(bool)), this, SLOT(togglePlot(bool)));
	QObject::connect(plotTRTSRadio, SIGNAL(toggled(bool)), splot, SLOT(setVisible(bool)));
	QObject::connect(plotTRTSRadio, SIGNAL(toggled(bool)), plotBox, SLOT(setVisible(bool)));
	QObject::connect(plotTRTSRadio, SIGNAL(toggled(bool)), this, SLOT(togglePlot(bool)));
	plotNoneRadio->setChecked(true);

	plotNoneRadio->setToolTip("Plot nothing");
	plotPRCRadio->setToolTip("Plot F1, F2");
	plotTRTSRadio->setToolTip("Plot tr-ts");
	
	QHBoxLayout *statusRowLayout = new QHBoxLayout;
	moduleStatus = new QLabel;
	moduleStatus->setFrameStyle(QFrame::Panel | QFrame::Sunken);
	moduleStatus->setAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
	moduleStatus->setText("Waiting...");
	statusRowLayout->addWidget(moduleStatus);

	QObject::connect(DefaultGUIModel::pauseButton, SIGNAL(toggled(bool)), savePlotButton, SLOT(setEnabled(bool)));
	QObject::connect(DefaultGUIModel::pauseButton, SIGNAL(toggled(bool)), saveDataButton, SLOT(setEnabled(bool)));
	QObject::connect(DefaultGUIModel::pauseButton, SIGNAL(toggled(bool)), printButton, SLOT(setEnabled(bool)));
	QObject::connect(DefaultGUIModel::pauseButton, SIGNAL(toggled(bool)), startBttn, SLOT(setEnabled(bool)));
	QObject::connect(DefaultGUIModel::pauseButton, SIGNAL(toggled(bool)), getPRCBttn, SLOT(setEnabled(bool)));
	QObject::connect(DefaultGUIModel::pauseButton, SIGNAL(toggled(bool)), DefaultGUIModel::modifyButton, SLOT(setEnabled(bool)));
	
	DefaultGUIModel::pauseButton->setToolTip("Start/Stop PRC protocol");
	DefaultGUIModel::modifyButton->setToolTip("Commit changes to parameter values");
	DefaultGUIModel::unloadButton->setToolTip("Close plug-in");
	QObject::connect(this, SIGNAL(PRC::setPlotRange(double, double, double, double)), splot, SLOT(setAxes(double, double, double, double)));

	customLayout->addWidget(baselineBox, 0, 0, 1, 1);
	customLayout->addWidget(plotBox, 0, 1, 1, 1);
	customLayout->addLayout(optionRowLayout, 2, 0, 1, 1);
	customLayout->addWidget(plotTypeBox, 3, 0, 1, 1);
	customLayout->addLayout(statusRowLayout, 4, 0, 1, 1);
	customLayout->addWidget(splot, 1, 1, 3, 1);

	plotBox->hide();
	splot->hide();

	setLayout(customLayout);
}

// custom functions
void PRC::initParameters() {
	dt = RT::System::getInstance()->getPeriod() * 1e-9; // period in s
	mindelay = 100e-3; // s
	maxdelay = 350e-3;
	stepsize = 100e-3;
	gmax = .01e-9; // S
	tau = 10e-3; // s
	esyn = 0; // V excitatory
	//  esyn = -70e-3;     // V inhibitory
	repeat = 1;
	randomize = false;
	plotPRC = false;
	saveData = false;
	
	// spike detector parameters
	thresh = -.02;
	min_int = 5e-3;
	spikestate = 0;
	last_spike = 0;
	
	bookkeep();
	srand(time(NULL));
	nstep = int((maxdelay + 1e-12 - mindelay) / stepsize) + 1; // calculate the number of steps
}

void PRC::initDelayArray() {
	delays = new double[nstep];
	
	for (int i = 0; i < nstep; i++) {
		delays[i] = mindelay + i * stepsize;
		delays[i] = delays[i] * 1000; // convert to ms
	}
	random = false;
	arrtr.assign(nstep, 0);
	arrts.assign(nstep, 0);
	arrP0.assign(nstep, 0);
	arrphi.assign(nstep, 0);
	arrPRC1.assign(nstep, 0);
	arrPRC2.assign(nstep, 0);
}

void PRC::shuffleonce() {
	if (stepcount == 0 && randomize == true && random == false) { // shuffle once per cycle
		std::random_shuffle(delays, delays + nstep);
		random = true;
	}
}

double PRC::shiftstimulus(double spiketime) {
	if (systime - spiketime >= delays[stepcount] / 1000) { // compute in s
		if (Gsyncount < tau * 10 / dt + 1) {
			Gsyncount++;
			return (Vm - esyn) * arrGsyn[Gsyncount - 1];
		} else {
			return 0;
		}
	} else {
		return 0;
	}
}

void PRC::initStimulus() {
	arrGsyn = new double[int(tau * 10 / dt + 1)]; // compute for 10 times tau to get long enough perturbation
	for (int i = 0; i < int(tau * 10 / dt + 1); i++) {
		arrGsyn[i] = -1 * gmax * (i * dt) / tau * exp(-(i * dt - tau) / tau);
	}
}

void PRC::bookkeep() {
	nidx = 0;
	count = 0; // reset time to zero
	stepcount = 0; // restart protocol
	cyclecount = 0;
	spikecount = 0;
	periodnum = 0;
	spikestate = 0;
	Gsyncount = 0;
	spktime = 0;
	prevspktime = 0;
	ISI = 0;
	P0 = 0;
	P1 = 0;
	P2 = 0;
	systime = 0;
	arrISI = new double[5];
	runningPeriod.clear();
}

void PRC::countspikes() {
	prevspktime = spktime;
	spktime = systime;
	periodnum = spikecount;
	spikecount++;
	ISI = (spktime - prevspktime) * 1000;
}

void PRC::resetcycle() {
	cyclecount++;
	stepcount = 0;
	random = false;
}

void PRC::makeplot() {
	if (plotGroup->checkedId() == 1) {
		if (arrphi.empty()) {
			printf("Error: nothing to plot\n");
		} else {
			for (int i = 0; i < (int) arrphi.size(); i++) {
				if (fabs(arrphi[i]) < 1e-12) arrphi[i] = 0;
				if (fabs(arrPRC1[i]) < 1e-12) arrPRC1[i] = 0;
				if (fabs(arrPRC2[i]) < 1e-12) arrPRC2[i] = 0;
			}
			emit PRC::setPlotRange(-.1, 1.1, -1.1, 1.1);
			emit drawData(&arrphi[0], &arrPRC1[0], (int) arrphi.size());
		}
	} else if (plotGroup->checkedId() == 2) {
		if (arrts.empty()) {
			printf("Error: nothing to plot\n");
		} else {
			double left = arrts[0];
			double right = arrts[0];
			double bottom = arrtr[0];
			double top = arrtr[0];
			for (int i = 1; i < (int) arrts.size(); i++) {
				left = std::min(left, arrts[i]);
				right = std::max(right, arrts[i]);
				bottom = std::min(bottom, arrtr[i]);
				top = std::max(top, arrtr[i]);
			}
			emit PRC::setPlotRange(left * 0.9, right * 1.1, bottom * 0.9, top * 1.1);
			//drawRefLines();
			emit drawData(&arrts[0], &arrtr[0], (int) arrts.size());
		}
	}
}

void PRC::clearData() {
	arrtr.clear(); // milliseconds
	arrts.clear(); // holds delays in order for each cycle in seconds
	arrP0.clear(); // milliseconds
	arrphi.clear(); // holds phases, unitless
	arrPRC1.clear(); // store 1st order PRC, unitless
	arrPRC2.clear(); // store 2nd order PRC, unitless
	printf("Data cleared.\n");
}

void PRC::print() {
#if 1
	QPrinter printer;
#else
	QPrinter printer(QPrinter::HighResolution);
#if QT_VERSION < 0x040000
	printer.setOutputToFile(true);
	printer.setOutputFileName("/tmp/PRC.ps");
	printer.setColorMode(QPrinter::Color);
#else
	printer.setOutputFileName("/tmp/PRC.pdf");
#endif
#endif
	QString docName = splot->title().text();
	if (!docName.isEmpty()) {
		docName.replace(QRegExp(QString::fromLatin1("\n")), tr(" -- "));
		printer.setDocName(docName);
	}
	printer.setCreator("PRC Curve");
	printer.setOrientation(QPrinter::Landscape);
#if QT_VERSION >= 0x040000
	QPrintDialog dialog(&printer);
	if (dialog.exec()) {
#else
	if (printer.setup()) {
#endif
/*
		RTXIPrintFilter filter;
		if (printer.colorMode() == QPrinter::GrayScale) {
			int options = QwtPlotPrintFilter::PrintAll;
			filter.setOptions(options);
			filter.color(QColor(29, 100, 141),
			QwtPlotPrintFilter::CanvasBackground);
			filter.color(Qt::white, QwtPlotPrintFilter::CurveSymbol);
		}
*/
//		splot->print(printer, filter);
		QwtPlotRenderer *renderer = new QwtPlotRenderer;
		renderer->renderTo(splot, printer);
	}
}
	
void PRC::exportSVG() {
	QString fileName = "PRC.svg";
#if QT_VERSION < 0x040000
#ifndef QT_NO_FILEDIALOG
	fileName = QFileDialog::getSaveFileName("PRC.svg", "SVG Documents (*.svg)", this);
#endif
	if (!fileName.isEmpty()) {
		// enable workaround for Qt3 misalignments
		QwtPainter::setSVGMode(true);
		QPicture picture;
		QPainter p(&picture);
		splot->print(&p, QRect(0, 0, 800, 600));
		p.end();
		picture.save(fileName, "svg");
	}
#elif QT_VERSION >= 0x040300
#ifdef QT_SVG_LIB
#ifndef QT_NO_FILEDIALOG
	fileName = QFileDialog::getSaveFileName(this, "Export File Name", QString(),"SVG Documents (*.svg)");
#endif
	if (!fileName.isEmpty()) {
		QSvgGenerator generator;
		generator.setFileName(fileName);
		generator.setSize(QSize(800, 600));
		splot->print(generator);
	}
#endif
#endif
}

void PRC::savePRC() {
	QFileDialog* fd = new QFileDialog(this, "Save File As"/*, TRUE*/);
	fd->setFileMode(QFileDialog::AnyFile);
	fd->setViewMode(QFileDialog::Detail);
	QString fileName;
	if (fd->exec() == QDialog::Accepted) {
		QStringList files = fd->selectedFiles();
		if (!files.isEmpty()) fileName = files.takeFirst();
//		fileName = fd->selectedFile();
		if (OpenFile(fileName)) {
			stream << "phase mean(P0) PRC1 PRC2 ts tr\n";
			for (int i = 0; i < (int) arrts.size(); i++) {
				stream << (double) arrphi[i] << " " << (double) arrP0[i] << " "
				<< (double) arrPRC1[i] << " " << (double) arrPRC2[i] << " "
				<< (double) arrts[i] << " " << (double) arrtr[i] << "\n";
			}
			dataFile.close();
			printf("File closed.\n");
		} else {
			QMessageBox::information(this, "PRC: Save PRC data",
			"There was an error writing to this file. You can view\n"
			"the computed values in the terminal.\n");
		}
	}
}

void PRC::startBaseline() {
	update(STARTBASELINE);
}

void PRC::measurePRC() {
	update(MEASUREPRC);
}

void PRC::togglerandom(bool on) {
	randomize = on;
}

void PRC::togglePlot(bool on) {
	plotPRC = on;
	if (on) makeplot();
	QTimer::singleShot(0, this, SLOT(resizeMe()));
}

bool PRC::OpenFile(QString FName) {
	dataFile.setFileName(FName);
	if (dataFile.exists()) {
		switch (QMessageBox::warning(this, "PRC", tr("This file already exists: %1.\n").arg(FName), 
		                             "Overwrite", "Append", "Cancel", 0, 2)) {

			case 0: // overwrite
				dataFile.remove();
				if (!dataFile.open(QIODevice::Unbuffered | QIODevice::WriteOnly)) return false; 
				break;

			case 1: // append
				if (!dataFile.open(QIODevice::Unbuffered | QIODevice::WriteOnly | QIODevice::Append)) return false;
				break;

			case 2: // cancel
				return false;
				break;
		}
	} else {
		if (!dataFile.open(QIODevice::Unbuffered | QIODevice::WriteOnly)) return false;
	}
	stream.setDevice(&dataFile);
	printf("File opened: %s\n", FName.toStdString().data());
	return true;
}
