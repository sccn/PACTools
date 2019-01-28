# EEGLAB Event Related PAC Tool
The Event Related PAC Tool (ERPAC Tool) is an EEGLAB plug-in to compute phase-amplitude coupling. 
While the plug-in is still in its early development, usage and feedback is encouraged.
In addition to traditional methods to compute PAC, the plugin include the Instantaneuous and Event-Related implementation of the Mutual Information Phase-Amplitude Coupling Method (MIPAC) (see Martinez-Cancino et al 2019).
The toolbox is developed and maintained at the Swartz Center for Computational Neurosciences, UCSD, La Jolla, California.


##Table of Contents
1. [Phase-amplitude coupling in neurosciences](#phase-amplitude-coupling-in-neurosciences)
2. [Methods implemented in the toolbox](#methods-implemented-in-the-toolbox)
3. [Plugin architecture and workflow](#plugin-architecture-and-workflow)
	1. [Plugin architecture]([#plugin-architecture])
	2. [Graphical user interface](#graphical-user-interface)
	3. [Output structure](#Output structure)
4. [Setting up the plug-in](#setting-up-the-plug-in)  
5. [Demos](#demos)
   1. [Demo 1: Using the GUI](#demo-1:-Using-the-GUI)
   2. [Demo 2: Mean Vector Length modulation Index in a single continuous signal](#demo-2:-mean-vectorlength-modulation-index-in-a-single-continuous-signal)
   3. [Demo 3: Instantaneous MIPAC in a single continuous signal]([#demo-3:-instantaneous-mipac-in-a-single-continuous-signal])
   4. [Demo 4: General Linear Model PAC in signals with multiple trials](#demo-4:-general-linear-model-pac-in-signals-with-multiple-trials)
   5. [Demo 5: Event-Related  MIPAC in signals with multiple trials](#demo-5:-event-related-mipac-in-signals-with-multiple-trials)
5. [Contributions and feedback](#contributions-and-feedback)

## Phase Amplitude Coupling in Neurosciences 
Cross-frequency coupling (CFC) could refer to any possible interaction between frequencies, phases and amplitudes of oscillatory phenomena (*Sotero, 2016*), most experimental work has focused on three types of CFC: amplitude-amplitude coupling (AAC) or comodulation, phase-phase coupling (PPC) including bicoherence, and phase-amplitude coupling (PAC). Among them, PAC has attracted increasing interest given the growing amount of evidence of its potential role in brain information processing and its changes under pathological conditions including epilepsy (*López-Azcárate et al., 2010; De Hemptinne et al., 2013*). In PAC, the instantaneous amplitude of a higher frequency band within a signal is modulated by (or otherwise linked to) the instantaneous phase of a lower-frequency band of the same (or a different) signal.



## Methods Implemented in the Toolbox
Several methods have been proposed to measure PAC. However, none is currently a gold standard.  In this toolbox, in addition to the recently developed Mutual Information Phase Amplitude Coupling (MIPAC) (*Martinez-Cancino et al., 2019*) , we have currently implemented three of the measures most often cited in the PAC literature: the Mean Vector Length Modulation Index (MVLmi) (*Canolty et al., 2006*), the Kullback-Leibler Modulation Index (KLmi) (*Tort et al., 2010*), and the General Linear Model Modulation Index (GLMmi) (*Penny et al., 2008a*). These measures have the ability to operate either in continues and epoched signals. In the case of epoched signals, a scheme similar to the one proposed by *Voytek et al., 2013* with the use of the method by (*Penny et al., 2008*) in the dimension of the trials (assuming a data matrix of dimensions equal to number of trials by latencies) is implemented.
 

### Continuous signal(s)
In the table below are listed the methods implemented to compute PAC in continuous signals. References to the specific methods are listed in the second rowof the table. The dimension of the output of the PAC methods  is indiceted in the third row of the table. 


Most of the PAC methods in the literature return a single value of the PAC measure for a given single trial signal (two single trials signals corresponding to different frequency bands of interest can be assumed as well). In the case of Instantaneous MIPAC, though, the result is provided as a unidimensional time series describing the PAC dynamics in the input signal(s). 

| Method                                | Reference                         | Output Dimension | Notes
| ---------                             | -----------                       | --------------   | ----
| Mean Vector Length Modulation Index   | [Canolty et al., 2006]()          | Single value     |
| Kullback-Leibler Modulation Index     | [Tort et al., 2010]()             | Single value     |
| General Linear Model Modulation Index | [Penny et al., 2008]()            | Single value     |
| Instantaneous Mutual Information PAC  | [Martinez-Cancino et al., 2019]() | Unidimensional   |

### Epoched signal(s)
In the table below are listed the current methods implemented in the toolbox to estimate PAC in epoched data. Epoched data is usually the result of extracting snipets of signals time-locked to an event(s) of interest. Here epoched data is assumed as being formated as a data matrix with dimensions of number of epochs(trials) by number of latencies(timepoints).  The three first methods listed in the table are a natural extension of the methods listed in the previous seccion but applying them onto each latency along the dimension of the epochs. The fisrt application of this scheme was proposed by *Voytek et al., 2013* as an extension of the method by Penny et al., 2008. These methods return a PAC time series describing the 'average' dynamics of the procces in the trials. Event-related MIPAC method, though, return a PAC time series for each trial provided.

| Method                                | Reference                         | Output Dimension | Notes
| ---------                             | -----------                       | --------------   | ----
| Mean Vector Length Modulation Index   | [Canolty et al., 2006]()          | Unidimensional   |
| Kullback-Leibler Modulation Index     | [Tort et al., 2010]()             | Unidimensional   |
| General Linear Model Modulation Index | [Voytek et al., 2013]()           | Unidimensional   |
| Event-Related Mutual Information PAC  | [Martinez-Cancino et al., 2019]() | Bidimensional    |


## Plugin Architecture and Workflow
### Plugin architecture
The plugin ERPAC is developed as an EEGLAB plugin. Given this, it shares the same philosophy regarding the functions structure and hierarchy as well as data formats(.set) as EEGLAB. Functions in EEGLAB are designed to provide users, both novice and expert Matlab users, with an easy and flexible usage. Depending on their level of Matlab expertise, users can either interact only with the graphics interface (GUI), else they can call functions directly from the Matlab command line or write their own Matlab scripts using EEGLAB functions and structures. This arrangement defines the hierarchy implemented by the two-level functions used in ERPAC toolbox. 
 
 Specifically in the plugin, the top layer function pop_pac.m provides its own GUI. Called with no (or few) arguments (as from the EEGLAB GUI), this function pops up a query window to gather additional parameter choices. The pop\_pac.m function can also be called directly from the Matlab command line or from Matlab scripts. 
 
As the top-layer function, pop_pac.m  provides the front-end interface for the toolbox it also serves as the bridge to the inner layer function, eeg\_pac.m. Users with a high level of Matlab/EEGLAB expertise can call this function directly by providing all the inputs required. The function eeg\_pac.m is indeed, the core function of the toolbox and is responsible for processing and parsing the input data and options in order to distribute it to the functions in charge of the computation of each of the PAC methods mentioned in the section [Methods Implemented in the Toolbox](#methods-implemented-in-the-toolbox).

### Graphical user interface
ERPAC tool provideds a flexible GUI that allow users to take advantage of the same functionalities provided from the command line. To invoke the GUI from the EEGLAB  GUI, click the menu *Tools >  ERPAC Tool > Estim. PAC*, otherwise you can launch the gui from the command line by typing `EEG = pop_pac(EEG);`. The figure below shows the graphical user interface of the toolbox. 
  
<Figure: GUI ===================================================================== >
<center>
<img  style="float: center;" src="doc/img/fig_gui_erpac_sample.tiff" alt="drawing" width="500"/>
<end>
</center>

The GUI is divided in four parts designated by the labels: **Data type**/**CFC type**, **PAC method**, **Command line options** and **PAC statistics**.
In the first section (**Data type**/**CFC type**), the type of data used for PAC computation can be selected in **Data type** between channel data (*Channels*) or ICA decomposed data (*Components*). The label **CFC type** is a static text indicating the type of CFC computation performed. The rationale of keeping this in the GUI is that future releases of the toolbox may contain other CFC methods in addition to PAC.
Right in the next line, a set of edits are used to input the property values for **Phase data**  and **Amplitude data**. In the first column, the index of the channels/components to use to compute PAC are defined (**Comp/chan indices**). In the second column (**Freq range [lo hi] Hz**), the range of frequencies (in Hz) to compute the instantaneous phase and amplitude can be defined. The number of frequencies in these ranges can be defined in the last column (**# Frequencies**.

The nex two sections allow for the selection of the PAC method (**PAC Method**) and input of optional parameters at **Command line options**. 

The last section comprises the settings for the computation of PAC statistics (**PAC statistics**). Here the number of surrogates (**# surrogates**), number of blocks to use to shuffle the data for generating the surrogates (**# blocks**), the significance threshold (**Significance threshold (0<p<1)**) and multiple comparison correction (**Correct for multiple comparisons**) can bet set. Three buttons lay at the bottom  of the GUI designated to launch the help documention (button: **Help**), cancel the execution of the GUI without further action (button: **Cancel**) and to start the execution of PAC computation with the settings provided (button: **OK**).
 
 
### Structure of outputs
 When computing PAC from pop\_pac.m,  the results of the computation are stored in the field *EEG.etc.pac.eegpac* structure. As an example, in the snippet below is shown a sample structure storing the computation of PAC using Instantaneous MIPAC (*instmipac*) and Kullback-Leibler Modulation Index methods.
 
```matlab
>> EEG.etc.eegpac

  struct with fields:

     chanindx: {[1 1]}
     chantype: 1
       params: [1×1 struct]
    instmipac: {[1×1 struct]}
         klmi: {[1×1 struct]}
```


The last two fields here indicate that the measure computed was Instantaneous MIPAC (*instmipac*). This field takes the name of the measure computed. In practice, we may find as many fields like this as PAC measures computed. This fields, in general, store the PAC values, the dimension of the output and results specific to the method computed.
 
### Signal Processing Prior to Computing PAC

 
## Demos
### Demo 1: Mean Vector Length Modulation Index and Instantaneous MIPAC in a single continuous signal

#### Computation

In this demo we will show how to compute PAC from the GUI. For this, lets first load the sample dataset *strial\_sim\_pac\_ph8amp60.set* into EEGLAB. The dataset cab be loaded either from the EEGLAB GUI or from the command windows by using the following code: 

```matlab
EEG = pop_loadset('filename','strial_sim_pac_ph8amp60.set');
eeglab redraw; 
```

Note: The current directory is assumed to be the folder containing the toolbox.

This dataset contains a simulated PAC signal where the instantaneous phase at 8Hz and the instantaneous amplitude at 60Hz are coupled during two segments of the signal, (see figure below).

<center>
<img src="doc/img/fig_sim_signal_1trial.tiff" alt="drawing" width="600"/>
</center>

After loading the dataset, we will proceed to compute PAC using ERPAC  from its main GUI. To launch the GUI, type `EEG = pop_pac(EEG)` in the MATLAB command windows (see below) and enter the parameters as shown in the figure below:

#### Visualization

`h = eeg_plotpac(EEG,'Comodulogram', 'pacmethod','mvlmi'`

`h = eeg_plotpac(EEG,'PhaseAmpTime', 'pacmethod','instmipac'`


The
### Demo 2: Computing Mean Vector Length Modulation Index and Event Related MIPAC in a signal with multiple trials
#### Computation

#### Visualization

### Contributions and feedback
This is an open source project, however, since this is still a Beta Version, please, contact the author for contributions.

