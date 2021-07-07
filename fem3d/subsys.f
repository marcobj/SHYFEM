
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-2006,2008-2020  Georg Umgiesser
!    Copyright (C) 2006-2008,2014-2015,2019  Christian Ferrarin
!    Copyright (C) 2008  Andrea Cucco
!    Copyright (C) 2019  Petras Zemlys
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

c routines system dependent
c
c contents :
c
c subroutine nlsinh             initializes the hp parameter file
c subroutine nlsina             initializes the ap parameter file
c subroutine fnminh             initializes default names
c
c revision log :
c
c 22.01.1998	ggu	no more strdir, apndir -> read apnpar.str from local
c 18.03.1998	ggu	no '[not defined...]' in basdir...
c 20.03.1998	ggu	for dos systems .memory -> _memory
c 22.03.1998	ggu	some changes for Technital integrated
c 30.04.1998	ggu	some changes for flux sections
c 13.05.1998	ggu	added colfil (apncol.str)
c 25.05.1998	ggu	documentation added
c 26.05.1998	ggu	documentation for wind file added
c 19.06.1998	ggu	some useless comments deleted
c 24.06.1998	ggu	qflux (heat flux) added
c 22.07.1998	ggu	more on documentation
c 23.07.1998	ggu	documentation
c 12.08.1998	ggu	new parameter dlat -> specify latitude for coriolis
c 02.09.1998	ggu	hack: change depth in Venice inlets (hlido,...)
c 24.11.1998	ggu	switch for biological reactor introduced
c 12.02.1999	ggu	default of parameter file changed to apnstd.str
c 12.02.1999	ggu	new parameters for plotting
c 13.04.1999	ggu	new parameter itide
c 19.04.1999	ggu	itide changed to rtide
c 27.05.1999	ggu	icust introduced
c 31.05.1999	ggu	dval default changed
c 28.09.1999	ggu	new flag regflx
c 19.11.1999	ggu	new parameters for section vol
c 08.08.2000	ggu	hlvmin is now percentage of last layer thickness
c 03.12.2001	ggu	new parameters (rlhdif)
c 11.10.2002	ggu	rstot new meaning
c 14.10.2002	ggu	atpar,adpar,aapar -> default set to 1.0 (implicit)
c 05.10.2003	ggu	changes in color handling of post routines
c 04.12.2003	ggu	sediment and wave module integrated
c 05.03.2004	ggu	bio variable for initialization
c 03.09.2004	ggu	restart variables
c 22.09.2004	ggu	nlsina_3d(): always use 3d file for plot (level=0)
c 28.09.2004	ggu	lagrangian routines integrated (LAGR)
c 05.10.2004	ggu	some more documentation
c 02.12.2004	ggu	documentation for variable time step
c 06.12.2004	ggu	new subroutine nlsina_legvar
c 17.01.2005	ggu	new parameters for horizontal diffusion
c 15.03.2005	ggu	cleaned horiz. diff., read new section legvar
c 19.05.2005	ggu	added time for custom routine (tcust)
c 04.11.2005	ggu	new parameter itlin (semi-lagrangian), some corrections
c 07.11.2005	ggu	new parameter itvd (TVD)
c 07.11.2005	ggu	parameters deleted: isedi,sedref,sedgrs
c 08.11.2005	ggu	do not set old parameters, some in nlsina_para
c 08.11.2005	ggu	more documentation
c 16.02.2006	ggu	new flag itoxi and file toxi
c 23.03.2006	ggu	new variable itunit for unit of time step
c 12.09.2006	ggu	time and date introduced
c 28.09.2006	ggu	new iprogr for progress of simulation
c 18.10.2006	ccf	new params iwave and dtwave for wave model
c 15.11.2006	ggu	new parameters to construct stretched vert. coordinates
c 16.11.2006	ggu	use itemp,isalt to decide about advection
c 29.11.2006	ggu	rwhpar for horizontal diffusion in lagrangian model
c 27.08.2007	ccf	isphe = 1  for spherical coordinate system
c 07.04.2008	aac	parameters for ersem ecological model call
c 10.04.2008	ccf	netcdf and gotmpa parameters
c 17.04.2008	ggu	new param ievap (to compute evaporation mass flux)
c 28.04.2008	ggu	rstot deleted, contau introduced
c 29.04.2008	ggu&aac	new vars for ERSEM
c 16.06.2008	ggu	old parts deleted, new documentation
c 17.09.2008	ggu	new interpretation for level (-1 = bottom)
c 11.10.2008	ggu	zfranco added
c 10.11.2008	ggu	new variable ilytyp, cleaned
c 19.11.2008	ggu	new variable noslip
c 24.11.2008	ggu	new variable vreps (mass error)
c 06.12.2008	ggu	new variables nbfix and nsigma
c 09.01.2009	ggu	documentation
c 21.01.2009	ggu	new variable vrerr (stop if mass error)
c 29.01.2009	ggu	various changes, better documentation
c 13.03.2009	ggu	bugfix: some parameters had default section not set
c 07.05.2009	ggu	new parameter ityrst
c 15.06.2009	ggu	new parameters for plot: isphe, reggrd, reggry
c 11.09.2009	ggu	new section $sect
c 09.10.2009	ggu	new parameter sclvel
c 13.10.2009	ggu	documentation is $sect, new hvmax, lvmax
c 22.02.2010	ggu	new parameter itdrag
c 23.02.2010	ggu	new parameter regdst
c 23.03.2010	ggu	changed v6.1.1
c 26.03.2010	ggu	new parameters for arrows in section plot
c 28.09.2010	ggu	new value for icor
c 29.09.2010	ggu	new param vmode,rxscal,ryscal
c 15.12.2010	ggu	nsigma renamed to nbsig, nsigma used for sigma layers
c 16.12.2010	ggu	changed VERS_6_1_15
c 21.12.2010	ggu	new parameter rwscal
c 27.01.2011	ggu	changed VERS_6_1_17
c 16.02.2011	ggu	new default for isphe, new routine fnm_aquabc_init()
c 25.02.2011	ggu	new param wsmax to catch errors in wind type
c 01.03.2011	ggu	changed VERS_6_1_20
c 23.03.2011	ggu	new parameter itvdv
c 24.03.2011	ggu	new parameters iheat,hdecay,botabs
c 14.04.2011	ggu	changed VERS_6_1_22
c 01.06.2011	ggu	new parameter idtmin
c 14.07.2011	ggu	changed VERS_6_1_27
c 18.08.2011	ggu	new parameter isoinp (interpolate inside element)
c 18.09.2011	ggu	change default for isphe for output (-1)
c 18.10.2011	ggu	changed VERS_6_1_33
c 03.11.2011	ggu	new parameter hsigma (hybrid)
c 18.11.2011	ggu	new subroutine nlsinh_proj() for projection
c 24.01.2012	ggu	new parameter nomp
c 14.02.2012	ggu	changed VERS_6_1_44
c 19.03.2012	ggu	changed VERS_6_1_49
c 02.05.2012	ggu	new default for ndccol (-> 0)
c 01.06.2012	ggu	changed VERS_6_1_53
c 29.08.2012	ggu	changed VERS_6_1_56
c 12.09.2012	ggu	changed VERS_6_1_57
c 24.10.2012	ggu	new parameter dxmin
c 10.05.2013	ggu	new parameters idtbox,itmbox, more comments
c 10.05.2013	ggu	new parameter inohyd
c 16.05.2013	ggu	file name bound renamed to zinit
c 13.06.2013	ggu	changed VERS_6_1_65
c 12.09.2013	ggu	changed VERS_6_1_67
c 25.10.2013	ggu	changed VERS_6_1_68
c 28.01.2014	ggu	changed VERS_6_1_71
c 07.03.2014	ggu	changed VERS_6_1_72
c 28.03.2014	ggu	some new params for lagrangian
c 10.04.2014	ccf	new section "wrt" for water renewal time
c 05.05.2014	ggu	changed VERS_6_1_74
c 30.05.2014	ggu	new default for dragco, new metpnt
c 18.07.2014	ggu	changed VERS_7_0_1
c 20.10.2014	ggu	new default for date (-1)
c 30.10.2014	ggu	changed VERS_7_0_4
c 05.11.2014	ggu	changed VERS_7_0_5
c 26.11.2014	ggu	changed VERS_7_0_7
c 01.12.2014	ccf	wave parameters moved to subwave.f
c 23.12.2014	ggu	changed VERS_7_0_11
c 09.01.2015	ggu	changed VERS_7_0_12
c 15.01.2015	ggu	changed VERS_7_1_1
c 26.02.2015	ggu	changed VERS_7_1_5
c 01.04.2015	ggu	changed VERS_7_1_7
c 30.04.2015	ggu	changed VERS_7_1_9
c 06.05.2015	ccf	new parameters itmoff and offlin
c 21.05.2015	ggu	changed VERS_7_1_11
c 29.09.2015	ggu	new boundary file surfvel
c 29.09.2015	ggu	new initial file uvinit, new flgrst
c 10.10.2015	ggu	changed VERS_7_3_2
c 19.10.2015	ggu	changed VERS_7_3_6
c 22.10.2015	ggu	changed VERS_7_3_8
c 05.11.2015	ggu	changed VERS_7_3_12
c 16.11.2015	ggu	changed VERS_7_3_14
c 20.11.2015	ggu	changed VERS_7_3_15
c 16.12.2015	ggu	changed VERS_7_3_16
c 01.02.2016	ggu	some plot params shifted to para section (bbgray, etc.)
c 19.02.2016	ggu	changed VERS_7_5_2
c 22.02.2016	ggu	new file for initial conditions bfmini
c 05.04.2016	ggu	new parameter iaicef for ice free areas
c 15.04.2016	ggu	changed VERS_7_5_8
c 25.05.2016	ggu	changed VERS_7_5_10
c 07.06.2016	ggu	changed VERS_7_5_12
c 10.06.2016	ggu	changed VERS_7_5_13
c 14.06.2016	ggu	changed VERS_7_5_14
c 17.06.2016	ggu	changed VERS_7_5_15
c 27.06.2016	ggu	changed VERS_7_5_16
c 09.09.2016	ggu	changed VERS_7_5_17
c 30.09.2016	ggu	changed VERS_7_5_18
c 05.10.2016	ggu	changed VERS_7_5_19
c 12.01.2017	ggu	changed VERS_7_5_21
c 13.02.2017	ggu	changed VERS_7_5_23
c 31.03.2017	ggu	changed VERS_7_5_24
c 13.04.2017	ggu	changed VERS_7_5_25
c 09.05.2017	ggu	changed VERS_7_5_26
c 26.05.2017	ggu	default of ishyff is 1 - no nos or ous files
c 13.06.2017	ggu	changed VERS_7_5_29
c 11.07.2017	ggu	changed VERS_7_5_30
c 26.09.2017	ggu	changed VERS_7_5_32
c 14.11.2017	ggu	changed VERS_7_5_36
c 05.12.2017	ggu	changed VERS_7_5_39
c 07.12.2017	ggu	changed VERS_7_5_40
c 24.01.2018	ggu	changed VERS_7_5_41
c 22.02.2018	ggu	changed VERS_7_5_42
c 03.04.2018	ggu	changed VERS_7_5_43
c 19.04.2018	ggu	changed VERS_7_5_45
c 06.07.2018	ggu	changed VERS_7_5_48
c 16.10.2018	ggu	changed VERS_7_5_50
c 25.10.2018	ggu	changed VERS_7_5_51
c 18.12.2018	ggu	changed VERS_7_5_52
c 27.12.2018	ggu	changed VERS_7_5_54
c 12.02.2019	ccf	introduced ibstrs for computing the bottom stess
c 16.02.2019	ggu	changed VERS_7_5_60
c 13.03.2019	ggu	new variables for plotting basin
c 21.05.2019	ggu	changed VERS_7_5_62
c 04.07.2019	ggu	new description for WW3
c 27.09.2019	pzy	new variables for aquabc
c 22.10.2019	ggu	some update of documentation
c 16.02.2020	ggu	itunit not supported anymore
c 05.03.2020	ggu	documentation upgraded
c 09.03.2020	ggu	more documentation upgraded, spell check
c 27.03.2020	ggu	documentation on ibarcl and nudging
c 11.11.2020	ggu	new parameter idtice and icemod
c 30.03.2021	ggu	new parameters idtdbg,itmdbg
c 23.04.2021	clr	new paramters petsc_zcfg and amgx_zcfg
c
c************************************************************************

	subroutine nlsinh

c initializes the parameter file for the main FE model

	implicit none

	call nlsinh_general
	call nlsinh_lagrg
        call nlsinh_wrt
	call nlsinh_bfmsc
	call nlsinh_proj
	call nlsinh_undoc
	call nlsinh_georg
	call nlsinh_unused
	call nlsinh_waves
	call nlsinh_nonhydro
	call nlsinh_connect

	end

c************************************************************************

	subroutine nlsinh_general

	use para

	implicit none

c $para section

	call sctpar('para')		!sets default section

c DOCS	START	S_para_h
c
c This section |$para| defines the general behavior of the simulation,
c gives various constants of parameters and determines what
c output files are written. In the following the meaning of
c all possible parameters is given.
c 
c Note that the only compulsory parameters in this section are
c the ones that chose the duration of the simulation and the
c integration time step. All other parameters are optional.
c
c DOCS	COMPULS		Compulsory time parameters
c
c These parameters are compulsory parameters that define the
c period of the simulation. They must be present in all cases.
c
c |itanf|	Start time of simulation. (Default 0)
c |itend|	End time of simulation.
c |idt|		Time step of integration.

	call addpar('itanf',0.)
	call addpar('itend',0.)
	call addpar('idt',0.)

cc------------------------------------------------------------------------

c DOCS	OUTPUT		Output parameters
c
c The following parameters deal with the output frequency
c and start time to external files. The content of the various
c output files should be looked up in the appropriate section.
c
c The default for the time step of output of the files is 0 which
c means that no output file is written. If the time step of the
c output files is equal to the time step of the simulation then
c at every time step the output file is written. The default start time
c of the output is |itanf|, the start of the simulation.
c
c |idtout|, |itmout|	Time step and start time for writing to file OUT,
c			the file containing the general hydrodynamic results.

	call addpar('idtout',0.)
	call addpar('itmout',-1.)

c |idtext|, |itmext|	Time step and start time for writing to file EXT,
c			the file containing hydrodynamic data of extra points.
c			The extra points for which the data is written
c			to this file are given in section |extra| of
c			the parameter file.

	call addpar('idtext',0.)
	call addpar('itmext',-1.)

c |idtrst|		Time step for writing the restart
c			file (extension RST). No restart file is written
c			with |idtrst| equal to 0. A negative value
c			is also possible for the time step. In this case
c			the time step used is |-idtrst|, but the file is
c			overwritten every time. It therefore contains 
c			always only the last written restart record. The
c			special value of |idtrst = -1| will write only the
c			last time step of the simulation in the restart file. 
c			This is useful if you want to start another
c			simulation from the last output. (Default 0)
c |itmrst|		Start time for writing the restart file. If
c			not given it is the beginning of the simulation.
c |itrst|		Time to use for the restart. If a restart
c			is performed, then the file name containing
c			the restart data has to be specified in |restrt|
c			and the time record corresponding to |itrst|
c			is used in this file. A value of -1 is also possible.
c			In this case the last record in the restart file
c			is used for the restart and the simulation starts
c			from this time. Be aware that a value of -1 changes
c			the parameter |itanf| to the time of the last
c			record found in |restrt|.
c |ityrst|		Type of restart. If 0 and the restart file is not
c			found the program will exit with an error. Otherwise
c			the program will simply continue with a cold start.
c			If |ityrst| is 1 and the given time record is not
c			found in the file it will exit with error. If
c			it is 2 it will initialize all values from the
c			first time record after |itrst|. Therefore, the
c			value of 2 will guarantee that the program will not
c			abort and continue running, but it might not
c			be doing what you intended. (Default 0)
c |flgrst|		This variable indicates which variables are read
c			from the restart file. By default all available
c			variables are read and used. If some variables
c			are not wanted (because, e.g., you want to restart
c			from a different T/S field), this fact can be indicated
c			in |flgrst|. 1 indicates restart of hydro values,
c			10 the depth values, 100 T/S values, 1000 the tracer
c			concentration, 10000 vertical velocities and
c			100000 the ecological variables. Therefore, a value
c			of 10111 indicates a restart of everything except
c			the tracer and the ecological values. The default
c			value for |flgrst| is -1, which means 111111.

	call addpar('idtrst',0.)
	call addpar('itmrst',-1.)
	call addpar('itrst',0.)
	call addpar('ityrst',0.)
	call addpar('flgrst',-1.)

c |idtres|, |itmres|	Time step and start time for writing to file RES,
c			the file containing residual hydrodynamic data.

	call addpar('idtres',0.)
	call addpar('itmres',-1.)

c |idtrms|, |itmrms|	Time step and start time for writing to file RMS,
c			the file containing hydrodynamic data of root mean
c			square velocities.

	call addpar('idtrms',0.)
	call addpar('itmrms',-1.)

c |idtflx|, |itmflx|	Time step and start time for writing to file FLX,
c			the file containing discharge data through defined
c			sections.
c			The transects for which the discharges are computed
c			are given in section |flux| of
c			the parameter file.

	call addpar('idtflx',0.)
	call addpar('itmflx',-1.)

c |idtstb|, |itmstb|	Time step and start time for writing a file with
c			the stability index for debug reasons.

	call addpar('idtstb',0.)
	call addpar('itmstb',-1.)

c |idtdbg|, |itmdbg|	Time step and start time for writing a debug file
c			with values of various variables. In order to
c			use this feature you have to run |shyfem| with
c			the command line option |-debout|. Please be
c			aware that this feature will create extremely
c			big files. Use at you own risk.

	call addpar('idtdbg',0.)
	call addpar('itmdbg',-1.)

c |idtvol|, |itmvol|	Time step and start time for writing to file VOL,
c			the file containing volume information of areas
c			defined by transects.
c			The transects that are used to compute the volumes
c			are given in section |volume| of
c			the parameter file.

	call addpar('idtvol',0.)
	call addpar('itmvol',-1.)

c |netcdf|		This parameter chooses output in NetCDF format 
c			if |netcdf| is 1, else the format is unformatted
c			FORTRAN files. (Default 0)

	call addpar('netcdf',0.)

c |idtoff|	handles offline mode (default 0):
c		\begin{description}
c		\item[0] do nothing (no offline routines called)
c		\item[$>$0] write offline data file (.off) 
c			with time step |idtoff|.
c		\item[$<$0] reads offline data from file |offlin|
c			defined in section |name|. Usage:
c			\begin{description}
c			\item[-1] uses offline hydro results
c			\item[-2] uses offline T/S results
c			\item[-4] uses offline turbulence results
c			\end{description}
c			Combinations are possible. A value of -3 would
c			read hydro and T/S results.
c		\end{description}

	call addpar('idtoff',0.)

c |itmoff|	Start time for writing to file OFF,
c		the file containing data for offline runs.
	call addpar('itmoff',-1.)

cc------------------------------------------------------------------------

c DOCS	TIME_DATE	General time and date parameters
c
c A time and date can be assigned to the simulation. These values
c refer to the time 0 of the FEM model. The format for the date is
c YYYYMMDD and for the time HHMMSS.
c You can also give a time zone if your time is not referring to
c GMT but to another time zone such as MET.

c |date|                The real date corresponding to time 0. (Default 0)
c |time|                The real time corresponding to time 0. (Default 0)
c |tz|                  The time zone you are in. This is 0 for GMT, 1 for MET
c                       and 2 for MEST (MET summer time). (Default 0)

        call addpar('date',-1.)
        call addpar('time',0.)
        call addpar('tz',0.)

cc------------------------------------------------------------------------

c DOCS	TERMS		Model parameters
c
c The next parameters define the inclusion or exclusion of
c certain terms of the primitive equations.
c
c |ilin|	Linearization of the momentum equations. If |ilin| 
c		is different from 0 the advective terms are not 
c		included in the computation. (Default 1)
c |itlin|	This parameter decides how the advective (non-linear)
c		terms are computed. The value of 0 (default) uses
c		the usual finite element discretization over a single
c		element. The value of 1 chooses a semi-lagrangian
c		approach that is theoretically stable also for
c		Courant numbers higher than 1. It is however recommended
c		that the time step is limited using |itsplt| and
c		|coumax| described below. (Default 0)
c |iclin|	Linearization of the continuity equation. If |iclin|
c		is different from 0 the depth term in the continuity
c		equation is taken to be constant. (Default 0)

	call addpar('ilin',1.)
	call addpar('itlin',0.)
	call addpar('iclin',0.)

cc undocumented: rlin can lower effect of non linear terms

	call addpar('rlin',1.)

c The next parameters allow for a variable time step in the
c hydrodynamic computations. This is especially important for the
c non-linear model (|ilin=0|) because in this case the criterion
c for stability cannot be determined a priori and in any case the
c time integration will not be unconditionally stable.
c
c The variable time steps allows for longer basic time steps
c (here called macro time steps) which have to be set in |idt|.
c It then computes the optimal time step (here micro time step)
c in order to not exceed the given Courant number. However,
c the value for the macro time step will never be exceeded.
c
c Normally time steps are always given in full seconds. This is still
c true when specifying the macro time step |idt|. In older versions also
c the computed micro time steps also had to be full integer values.
c Starting from version 7.1 also fractional time steps are allowed.
c This gives the possibility to have time steps smaller than 1 s.
c
c |itsplt|	Type of variable time step computation. If this value
c		is 0, the time step will be kept constant at its initial
c		value. A value of 1 divides the initial time step into
c		(possibly) equal parts, but makes sure that at the end
c		of the micro time steps one complete macro time
c		step has been executed. The mode |itsplt| = 2
c		does not care about the macro time step, but always
c		uses the biggest time step possible. In this case
c		it is not assured that after some micro time steps
c		a macro time step will be recovered. Please note
c		that the initial macro time step will never be exceeded.
c		In any case, the time step will always be rounded to the
c		next lower integer value. This is not the case with
c		|itsplt| = 3 where the highest possible fractional time step 
c		will be used. (Default 0)
c |coumax|	Normally the time step is computed in order to not
c		exceed the Courant number of 1. However, in some cases
c		the non-linear terms are stable even for a value higher
c		than 1 or there is a need to achieve a lower Courant number.
c		Setting |coumax| to the desired Courant number
c		achieves exactly this effect. (Default 1)
c |idtsyn|	In case of |itsplt| = 2 this parameter makes sure that
c		after a time of |idtsyn| the time step will be synchronized
c		to this time. Therefore, setting |idtsyn| = 3600 means
c		that there will be a time stamp every hour, even if the model
c		has to take one very small time step in order to reach that
c		time. This parameter is useful
c		only for |itsplt| = 2 and its default value of
c		0 does not make any synchronization.
c |idtmin|	This variable defines the smallest time step possible
c		when time step splitting is enabled. Normally the smallest
c		time step is 1 second. Please set |idtmin| to values
c		smaller than 1 in order to allow for fractional time steps.
c		A value of 0.001 allows for time steps of down to
c		1 millisecond. (Default 1)

	call addpar('itsplt',0.)
	call addpar('coumax',1.)
	call addpar('idtsyn',0.)
	call addpar('idtmin',1.)
	call addpar('tfact',0.)		!still to comment FIXME

c These parameters define the weighting of time level in the 
c semi-implicit algorithm. With these parameters the damping
c of gravity or Rossby waves can be controlled. Only modify them if
c you know what you are doing.
c
c |azpar|	Weighting of the new time level of the transport
c		terms in the continuity equation. (Default 0.5)
c |ampar|	Weighting of the new time level of the pressure
c		term in the momentum equations. (Default 0.5)
c |afpar|	Weighting of the new time level of the Coriolis
c		term in the momentum equations. (Default 0.5)
c |avpar|	Weighting of the new time level of the non-linear
c		advective terms in the momentum equations. (Default 0.0)

	call addpar('azpar',0.5)
	call addpar('ampar',0.5)
	call addpar('afpar',0.5)
	call addpar('avpar',0.0)

c The next parameters define the weighting of time level for the
c vertical stress and advection terms. They guarantee the stability
c of the vertical system. For this reason they are normally set to
c 1 which corresponds to a fully implicit discretization. Only
c modify them if you know what you are doing.
c
c |atpar|	Weighting of the new time level of the vertical
c		viscosity in the momentum equation. (Default 1.0)
c |adpar|	Weighting of the new time level of the vertical
c		diffusion in the scalar equations. (Default 1.0)
c |aapar|	Weighting of the new time level of the vertical
c		advection in the scalar equations. (Default 1.0)

	call addpar('atpar',1.0)	!time weighting for vertical viscosity
	call addpar('adpar',1.0)	!time weighting for vertical diffusion
	call addpar('aapar',1.0)	!time weighting for vertical advection

cc------------------------------------------------------------------------

c DOCS	CORIOLIS	Coriolis parameters
c
c The next parameters define the parameters to be used
c with the Coriolis terms.
c
c |icor|	If this parameter is 0, the Coriolis terms are
c		not included in the computation. A value of 1
c		uses a beta-plane approximation with a variable
c		Coriolis parameter $f$, whereas a value of
c		2 uses an f-plane approximation where the
c		Coriolis parameter $f$ is kept constant over the
c		whole domain. (Default 0)
c |dlat|	Average latitude of the basin. This is used to
c		compute the Coriolis parameter $f$. This parameter 
c		is not used if spherical coordinates are used 
c		(|isphe|=1) or if a coordinate 	projection is set 
c		(|iproj| $>$0). (Default 0)
c |isphe|	If 0 a cartesian coordinate system is used,
c		if 1 the coordinates are in the spherical system (lat/lon).
c		Please note that in case of spherical coordinates the
c		Coriolis term is always included in the computation, even
c		with |icor| = 0. If you really do not want to use the
c		Coriolis term, then please set |icor| = -1. The default is
c		-1, which means that the type of coordinate system will 
c		be determined automatically.

	call addpar('icor',0.)
	call addpar('dlat',100.)
	call addpar('isphe',-1.)

cc------------------------------------------------------------------------

c DOCS	DEPTH		Depth parameters
c
c The next parameters deal with handling depth values of the basin.
c
c |href|	Reference depth. If the depth values of the basin and
c		the water levels are referred to mean sea level, |href|
c		should be 0 (default value). Else this value is
c		subtracted from the given depth values. For example,
c		if |href = 0.20| then a depth value in the basin
c		of 1 meter will be reduced to 80 centimeters.

	call addpar('href',0.)

c |hzmin|	Minimum total water depth that will remain in a
c		node if the element becomes dry. (Default 0.01 m)
c |hzoff|	Total water depth at which an element will be
c		taken out of the computation because it becomes dry.
c		(Default 0.05 m)
c |hzon|	Total water depth at which a dry element will be
c		re-inserted into the computation.
c		(Default 0.10 m)

	call addpar('hzmin',0.01)
	call addpar('hzoff',0.05)
	call addpar('hzon',0.10)

c |hmin|	Minimum water depth (most shallow) for the whole
c		basin. All depth values of the basin will be adjusted
c		so that no water depth is shallower than |hmin|.
c		(Default is no adjustment)
c |hmax|	Maximum water depth (deepest) for the whole
c		basin. All depth values of the basin will be adjusted
c		so that no water depth is deeper than |hmax|.
c		(Default is no adjustment)

	call addpar('hmin',-99999.)
	call addpar('hmax',+99999.)

cc------------------------------------------------------------------------
cc friction - description in file subn35.f
cc------------------------------------------------------------------------

c \input{P_friction.tex}

	call addpar('ireib',0.)
	call addpar('czdef',0.)
	call addpar('iczv',1.)
	call addpar('uvmin',0.2) !Minimum current speed for ireib=10

cc------------------------------------------------------------------------

c DOCS	PHYSICAL		Physical parameters
c
c The next parameters describe physical values that can be adjusted
c if needed.
c
c |rowass|	Average density of sea water. (Default 1025 \densityunit)
c |roluft|	Average density of air. (Default 1.225 \densityunit)
c |grav|	Gravitational acceleration. (Default 9.81 \accelunit)

        call addpar('rowass',1025.)
        call addpar('roluft',1.225)
        call addpar('grav',9.81)

cc------------------------------------------------------------------------
cc wind - description in file submeteo2.f
cc------------------------------------------------------------------------

c \input{P_wind.tex}

        call addpar('iwtype',1.)
	call addpar('itdrag',0.)
	call addpar('dragco',2.5e-3)
	call addpar('wsmax',50.)
	call addpar('wslim',-1.)

cc------------------------------------------------------------------------

c DOCS  meteo                      Meteo and heat flux parameters

c The next parameters deal with the heat and meteo forcing.

c |iheat|	The type of heat flux algorithm (Default 1):
c		\begin{description}
c		\item[1] As in the AREG model
c		\item[2] As in the POM model
c		\item[3] Following A. Gill
c		\item[4] Following Dejak
c		\item[5] As in the GOTM model
c		\item[6] Using the COARE3.0 module
c		\item[7] Read heat fluxes directly from file. The columns
c		in the data file must be "time srad qsens qlat qlong".
c		\item[8] MFS heat fluxes as Pettenuzzo et al., 2010
c		\end{description}
c		Except when |iheat| is 7, the time series file has
c		the columns "time srad airt rhum cc".

	call addpar('iheat',1.)		!type of heat flux routine

c |ihtype|	Different ways of how to specify water vapor content
c		are possible. Normally relative humidity has to be
c		given (|ihtype|=1). However, also wet bulb temperature
c		(|ihtype|=2), dew point temperature (|ihtype|=3), or
c		specific humidity (|ihtype|=4) can
c		be given. (Default 1).

	call addpar('ihtype',1.)	!type of water vapor

c |isolp|	The type of solar penetration parameterization
c		by one or more exponential decay curves.
c		|isolp| = 0 sets an e-folding decay of radiation 
c		(one exponential decay curve) as function of depth |hdecay|.
c		|isolp| = 1 sets a profile of solar radiation with two length
c		scale of penetration. Following the Jerlov (Jerlov, N. G., 1968
c		Optical Oceanography, Elsevier, 194pp) classification the type
c		of water is clear water (type I). 
c		(Default 0) 

        call addpar('isolp',0.)         !type of solar penetration   

c |iwtyp|	The water types from clear water (type I) to the most 
c		turbid water (coastal water 9) following the classification of
c		Jerlov (Jerlov, N. G., 1968 Optical Oceanography, 
c		Elsevier, 194pp). The possible values for |iwtyp| are:
c		\begin{description}
c		\item[0] clear water type I 
c		\item[1] type IA
c		\item[2] type IB
c		\item[3] type II
c		\item[4] type III
c		\item[5] type 1
c		\item[6] type 3
c		\item[7] type 5
c		\item[8] type 7
c		\item[9] type 9
c		\end{description}

        call addpar('iwtyp',0.)         !water type (works if |isolp|=1)

c |hdecay|	Depth of e-folding decay of radiation [m]. If |hdecay| = 0 
c		everything is absorbed in first layer (Default 0).

	call addpar('hdecay',0.)	!depth of e-folding decay of radiation

c |botabs|	Heat absorption at bottom [fraction] (Default 0).
c		\begin{description}
c		\item[=0] everything is absorbed in last layer
c		\item[=1] bottom absorbs remaining radiation
c		\end{description}

	call addpar('botabs',0.)	!heat absorption at bottom

c |albedo|	General albedo (Default 0.06).

	call addpar('albedo',0.06)	!general albedo

c |albed4|	Albedo for temp below 4 degrees (Default 0.06).

	call addpar('albed4',0.06)	!albedo for temp below 4 degrees

c |ievap| 	Compute evaporation mass flux (Default 0).

	call addpar('ievap',0.)		!compute evaporation mass flux

cc------------------------------------------------------------------------

c DOCS	3D			Parameters for 3d
c
c The next parameters deal with the layer structure in 3D.

c |dzreg|	Normally the bottom of the various layers are given in
c		section |\$levels|. If only a regular vertical grid is desired
c		then the parameter |dzreg| can be used. It specifies the spacing
c		of the vertical layers in meters. (Default is 0, which means 
c		that the layers are specified explicitly in |\$levels|.

	call addpar('dzreg',0.)		!regular vertical grid

c The last layer (bottom layer) is treated in a special way. Depending on
c the parameter |ilytyp| there are various cases to be considered. A value
c of 0 leaves the last layer as it is, even if the thickness is very small.
c A value of 1 will always eliminate the last layer, if it has not full
c layer thickness. A value of 2 will do the same, but only if the last layer
c is smaller than |hlvmin| (in units of fraction). Finally, a value of
c 3 will add the last layer to the layer above, if its layer thickness
c is smaller than |hlvmin|.
c
c |ilytyp|	Treatment of last (bottom) layer. 0 means no adjustment,
c		1 deletes the last layer, if it is not a full layer,
c		2 only deletes it
c		if the layer thickness is less than |hlvmin|, and 3
c		adds the layer thickness to the layer above if it is smaller
c		than |hlvmin|. Therefore, 1 and 2 might change the
c		total depth and layer structure, while 3 only might
c		change the layer structure. The value of 1 will always
c		give you full layers at the bottom. (Default 3)
c |hlvmin|	Minimum layer thickness for last (bottom) layer used when
c		|ilytyp| is 2 or 3. The unit is fractions of the nominal
c		layer thickness. Therefore, a value of 0.5 indicates that
c		the last layer should be at least half of the full
c		layer. (Default 0.25)

	call addpar('ilytyp',3.00)	!type of depth adjustment
	call addpar('hlvmin',0.25)	!min percentage of last layer thickness

c The above parameters are dealing with zeta layers, where every layer
c has constant thickness, except the surface layer which is varying with
c the water level. The next parameters deal with sigma layers where all
c layers have varying thickness with the water level.

c |nsigma|	Number of sigma layers for the run. This parameter can
c		be given in alternative to specifying the sigma layers
c		in |\$levels|. Only regularly spaced sigma levels
c		will be created. (Default 0)
c |hsigma|	This is still an experimental feature. It allows to use
c		sigma layers above zeta layers. |hsigma| is the depth where
c		the transition between these two types of layers
c		is occurring. (Default 10000)

	call addpar('nsigma',0.)	!number of sigma layers
	call addpar('hsigma',10000.)	!lower depth of sigma layers (hybrid)

c The next parameters deal with vertical diffusivity and viscosity.

c |difmol|	Vertical molecular diffusivity parameter for
c		temperature, salinity, and tracer.
c		(Default 1.0e-06)
c |diftur|	Vertical turbulent diffusivity parameter for 
c		temperature, salinity, and tracer.
c		(Default 0)
c |vismol|	Vertical molecular viscosity parameter for momentum.
c		(Default 1.0e-06)
c |vistur|	Vertical turbulent viscosity parameter for momentum.
c		(Default 0)

        call addpar('difmol',1.0e-06)	!molecular vertical diffusivity
	call addpar('diftur',0.)	!diffusion parameter (vertical), cvpar
        call addpar('vismol',1.0e-06)	!molecular vertical viscosity
        call addpar('vistur',0.)	!turbulent vertical viscosity (nau)

c Instead of setting fixed values for viscosity and diffusivity, these
c parameters can also be computed by a turbulence closure scheme. 
c The parameter |iturb| defines what turbulence scheme is going to be used.
c A value of 0 uses no turbulence scheme. In this case be sure that 
c |vistur| and |diftur| have been set manually. |iturb=1| uses the GOTM
c routines for turbulence. In order to use this value |SHYFEM| has to be
c compiled with GOTM support. |iturb=2| uses a k-epsilon module, and a value
c of 3 uses the Munk-Anderson model. The recommended value for |iturb| 
c is 1 (GOTM module). (Default 0)

	call addpar('iturb',0.)		!use turbulence closure scheme

cc------------------------------------------------------------------------

c The next parameters deal with horizontal diffusion.
cc horizontal diffusion (Smagorinsky)

	call addpar('idhtyp',0.)	!type of horizontal diffusion/viscosity
	call addpar('dhlen',1.) 	!length scale for idhtyp=1
	call addpar('noslip',0.) 	!no slip conditions on boundary

	call addpar('idtype',2.) 	!type of hor diffus (delete after test)

	call addpar('ahpar',0.)		!austausch coefficient (still same??)

c |dhpar|	Horizontal diffusion parameter (general).
c		(Default 0)

	call addpar('dhpar',0.)		!diffusion parameter

c The next parameters deal with the control of the scalar transport 
c and diffusion equation. You have possibility to prescribe the tvd scheme
c desired and to limit the Courant number.
c
c |itvd|	Type of the horizontal advection scheme used for 
c		the transport and diffusion
c		equation. Normally an upwind scheme is used (0), but setting
c		the parameter |itvd| to a value greater than 0 
c		choses a TVD scheme. A value of 1 will use a TVD scheme
c		based on the average gradient, and a value of 2 will use
c		the gradient of the upwind node (recommended).
c		This feature is still experimental, so use with care. 
c		(Default 0)
c |itvdv|	Type of the vertical advection scheme used for 
c		the transport and diffusion
c		equation. Normally an upwind scheme is used (0), but setting
c		the parameter |itvdv| to 1 choses a TVD scheme. This feature
c		is still experimental, so use with care. (Default 0)
c |rstol|	Normally the internal time step for scalar advection is
c		automatically adjusted to produce a Courant number of 1
c		(marginal stability). You can set |rstol| to a smaller value 
c		if you think there are stability problems. (Default 1)

	call addpar('itvd',0.)		!horizontal TVD scheme?
	call addpar('itvdv',0.)		!vertical TVD scheme?
	call addpar('rstol',1.)		!limit time step to this Courant number

cc------------------------------------------------------------------------

c DOCS	VAR		Various parameters
c
c The next parameters describe various parameters not related to
c the above parameters.

c |tauvel|	If you have velocity observations given in file
c		|surfvel| then you can specify the relaxation
c		parameter $\tau$ in the variable |tauvel|. (Default 0,
c		which means no assimilation of velocities)

	call addpar('tauvel',0.)	!horizontal TVD scheme?

c |rtide|	If |rtide| = 1 the model calculates equilibrium tidal 
c		potential and load tides and uses these to force the 
c		free surface (Default 0).

	call addpar('rtide',0.)

c |ltidec|	Calibration factor for calculating the loading tide, 
c		which is computed in function of the total water depth as
c		$\beta=ltidec*H$. Usually it has a value of order 1e-6. 
c		If 0 no loading tide is computed (Default 0).

	call addpar('ltidec',0.)

c |ibstrs|	Call parameter for the routine bstress. If equal to 1 it 
c		computes (in function of currents and waves) and writes 
c		the bottom shear stess into a .bstress.shy file.
c		If 0 no bottom stess is computed (Default 0).

	call addpar('ibstrs',0.)

cc------------------------------------------------------------------------

c DOCS	ST	Temperature and salinity
c
c The next parameters deal with the transport and diffusion
c of temperature and salinity (T/S). 

c |ibarcl|	In order to compute
c		T/S the parameter |ibarcl| must be different from 0. 
c		Different values indoicate different ways to compute 
c		the advection and diffusion of the variables T/S 
c		and how their values is used in
c		the momentum equation. (Default 0)
c		\begin{description}
c		\item[0] No computation of T/S.
c		\item[1] Computation of T/S. The T/S field has a feedback
c			 through the baroclinic term on the momentum
c			 equation. This corresponds to a full
c			 baroclinic model.
c		\item[2] The model runs in diagnostic mode. No advection
c			 of the T/S field is done, but the baroclinic term
c			 is computed and used in the momentum equations.
c			 For this to work T/S fields have to be provided
c			 in external files |tempobs| and |saltobs|.
c		\item[3] The model computes the advection and diffusion
c			 of T/S, but their value is not used in the baroclinic
c			 terms in momentum equation. This corresponds to
c			 treating T/S as tracers.
c		\item[4] The model runs in baroclinic mode, but
c			 uses nudging for a basic data assimilation of T/S.
c			 For this to work T/S fields have to be provided
c			 in external files |tempobs| and |saltobs|.
c			 Moreover, a time scale for the nudging has to be
c			 provided, either as a constant using
c			 |temptaup| and |salttaip| or in spatially and
c			 temporarily varying fields given in external
c			 files |temptau| and |salttau|.
c		\end{description}

	call addpar('ibarcl',0.)	!compute baroclinic contributions ?

c In case |ibarcl| is different from 0 the variables T/S will
c be computed. However, they may be selectively turned off setting
c one of the two parameters |itemp| or |isalt| explicitly to 0.
c
c |itemp|	Flag if the computation on the temperature is done.
c		A value different from 0 computes the transport
c		and diffusion of the temperature. (Default 1)
c |isalt|	Flag if the computation on the salinity is done.
c		A value different from 0 computes the transport
c		and diffusion of the salinity. (Default 1)

	call addpar('itemp',1.)		!compute temperature ?
	call addpar('isalt',1.)		!compute salinity ?

c The next parameters set the initial conditions for temperature and salinity.
c Both the average value and and a stratification can be specified.
c Initial conditions can also be given in external files |tempin| and |saltin|.
c
c |temref|	Reference (initial) temperature of the water in
c		centigrade. (Default 0)
c |salref|	Reference (initial) salinity of the water in
c		psu (practical salinity units) or ppt.
c		(Default 0)
c |tstrat|	Initial temperature stratification in units of [C/km].
c		A positive value indicates a stable stratification.
c		(Default 0)
c |sstrat|	Initial salinity stratification in units of [psu/km].
c		A positive value indicates a stable stratification.
c		(Default 0)

	call addpar('temref',0.)	!reference temperature for T/S runs
	call addpar('salref',0.)	!reference salinity for T/S runs

	call addpar('sstrat',0.)	!salt stratification
	call addpar('tstrat',0.)	!temp stratification

c The next parameters deal with horizontal diffusion of temperature
c and salinity. These parameters overwrite the general parameter for
c horizontal diffusion |dhpar|.

c |thpar|	Horizontal diffusion parameter for temperature.
c		(Default 0)
c |shpar|	Horizontal diffusion parameter for salinity.
c		(Default 0)

	call addpar('thpar',-1.)	!horiz. diff. coeff. for temp.
	call addpar('shpar',-1.)	!horiz. diff. coeff. for sal.

c When using nudging (|ibarcl=4|) time scales for nudging have to given.
c They can be given in the parameters |temptaup| and |salttaup| or
c in external files |temptau| and |salttau|.

c |temptaup|	Time scale for temperature nudging in seconds.
c		(Default 0)
c |salttaup|	Time scale for salinity nudging in seconds.
c		(Default 0)

	call addpar('temptaup',-1.)	!tau for temperature
	call addpar('salttaup',-1.)	!tau for salinity

cc------------------------------------------------------------------------

c DOCS	CC	Concentrations
c
c The next parameters deal with the transport and diffusion
c of a conservative substance. The substance is dissolved in
c the water and acts like a tracer.
c
c |iconz|	Flag if the computation on the tracer is done.
c		A value different from 0 computes the transport
c		and diffusion of the substance. If greater than 1
c		|iconz| concentrations are simulated. (Default 0)
cc |iconza|	The concentration are normally written as instantaneous
cc		values. If |iconza| is different from 0, instead of
cc		the instantaneous values the minimum, average and maximm
cc		values are written. (Default 0)
c |conref|	Reference (initial) concentration of the tracer in
c		any unit. (Default 0)
c |taupar|	Decay rate for concentration if different from 0. In
c		this case |taupar| is the decay rate (e-folding time) in days.
c		This parameter is also used for multi-concentration runs.
c		In this case either one value has to be given that is used
c		for all concentrations, or |iconz| values have to be given,
c		one for each concentration.
c		(Default 0)
c |idecay|	Type of decay used. If 0 no decay is used.
c		A value of 1 uses the value of |taupar| as exponential decay.
c		A value of 2 uses a formulation of Chapra, where the
c		decay rate depends on T,S,light and settling. In this
c		case the value of |taupar| is ignored.
c		(Default 0)

	call addpar('iconz',0.)		!compute concentration ?
cc	call addpar('iconza',0.)	!average values ? - not working
	call addpar('conref',0.)	!reference concentration
	call addpar('idecay',0.)	!type of decay

	call para_deprecate('contau','taupar')	!no contau anymore -> use taupar
	call para_add_array_value('taupar',0.)	!decay rate [days]

c |chpar|	Horizontal diffusion parameter for the tracer.
c		This value overwrites the general parameter for
c		horizontal diffusion |dhpar|. (Default 0)

	call addpar('chpar',-1.)	!diffusion parameter

cc------------------------------------------------------------------------

c DOCS	STOUTPUT	Output for scalars
c
c The next parameters define the output frequency of the
c computed scalars (temperature, salinity, generic concentration) to file.
c
c |idtcon|, |itmcon|	Time step and start time for writing to file
c			conz.shy (concentration) and ts.shy (temperature and
c			salinity).
c |irho|	Flag if the density is also written together with T/S.
c		A value different from 0 writes the density to file.
c		(Default 0)

	call addpar('idtcon',0.)	!time step for output
	call addpar('itmcon',-1.)	!minimum time for output
	call addpar('irho',0.)		!write rho ?

c DOCS	END

cc------------------------------------------------------------------------
cc------------------------------------------------------------------------
cc------------------------------------------------------------------------
cc------------------------------------------------------------------------
cc------------------------------------------------------------------------

cc------------------------------------------------------------------------
cc still to be commented below here
cc------------------------------------------------------------------------

	call addpar('idtsti',0.)	!time step for stability index
	call addpar('itmsti',-1.)	!minimum time for stability index

cc------------------------------------------------------------------------

	call addpar('ihydro',1.)	!compute hydrodynamics (for debug)

cc------------------------------------------------------------------------

	call addpar('ibio',0.)		!run biological reactor
        call addpar('imerc',0.)		!run mercury reactor
        call addpar('issedi',0.)	!run simple sediments
	call addpar('itoxi',0.)		!run toxicological routines

cc------------------------------------------------------------------------

cc call routines in flux.f

	call addpar('regflx',0.)	!compute regular fluxes

cc------------------------------------------------------------------------

cc custom call

	call addpar('icust',0.)		!call custom routine
	call addpar('tcust',0.)		!time for custom routine

	call addpar('ishyff',1.)	!SHYFEM file format
					!0=old 1=new 2=both

	call addpar('ipoiss',0.)	!solve poisson equation

cc rain

	call addpar('zdist',0.)		!distributed water level

	call addpar('idtbox',0.)	!for boxes
	call addpar('itmbox',-1.)

	call addpar('ihwadv',0.)        !vert advection of horiz momentum

cc ice

	call addpar('iaicef',-99.)	!area code for ice free condition
	call addpar('icemod',0.)	!use ice model?
	call addpar('idtice',0.)	!time step for ice model

	end

c************************************************************************

	subroutine nlsinh_lagrg

	implicit none

c $lagrg section

c DOCS	START	P_lagrg
c
c Section |$lagrg| describes the use of the Lagrangian Particle Module.
c
c The lagrangian particles can be released:
c \begin{itemize}
c \item inside the given areas (filename |lgrlin|). If this file is not 
c      specified they are released over the whole domain. The amount of
c      particles released and the time step are specified by |nbdy| and 
c      |idtl|.
c \item at selected times and location, e.g. along a drifter track
c      (filename |lgrtrj|). |nbdy| particles are released at the times
c      and location specified in the file.
c \item as initial particle distribution (filename |lgrini|) at time
c      |itlgin|. This file has the same format as the lagrangian output.
c \item at the open boundaries either as particles per second or per
c      volume flux (parameter |lgrpps|).
c \end{itemize}
c
c |lgrlin| and |lgrtrj| are mutually exclusive. 
c
c The lagrangian module runs between the times |itlanf| and |itlend|. If one
c or both are missing, the simulation extremes are substituted. Inside
c the lagrangian simulation window, the release of particles inside a given
c area is controlled by the parameters |idtl|, |itranf| and |itrend|. 
c |itranf| gives the time of the first release, |itrend| the time for 
c the last release. If not given they are set equal to the extremes of 
c the lagrangian simulation. |idtl| is giving the time step of release.
c
c The output frequency of the results can be controlled by 
c |idtlgr| and |itmlgr|.
c
c Please find all details here below.

        call sctpar('lagrg')             !sets default section
        call sctfnm('lagrg')

c |ilagr|	Switch that indicates that the lagrangian module
c		should be run (default 0):
c		\begin{description}
c		\item[0] do nothing
c		\item[1] surface lagrangian
c		\item[2] 2d lagrangian
c		\item[3] 3d lagrangian
c		\end{description}

	call addpar('ilagr',0.)         !LAGR

c |nbdymax|	Maximum numbers of particles that can be in the domain.
c		This should be the maximum number of particles
c		that can be created and inserted. Use 0 to not limit
c		the number of particles (on your own risk). This
c		parameter must be set and has no default.

        call addpar('nbdymax',-1.)

c |nbdy|	Total numbers of particles to be released in the domain each
c		time a release of particles takes place. 
c		(Default 0)

        call addpar('nbdy',0.)

c |rwhpar|	A horizontal diffusion can be defined for the lagrangian model.
c		Its value can be specified in |rwhpar| and the units are 
c		[m**2/s]. If |rwhpar| < 0 the diffusion parameter depends on
c		the local diffusivity (see |idhtyp|) (Default 0)

	call addpar('rwhpar',0.)	!diffusion for lagrangian model

c |itlanf, itlend|	The start and end time for the lagrangian module.
c			If not given, the module runs for the whole simulation.

        call addpar('itlanf',-1.)
        call addpar('itlend',-1.)

c |itmlgr, idtlgr|	Initial time and time step for the output to file
c			of the particles. if |idtlgr| is 0, 
c			no output is written. (Default 0)

        call addpar('itmlgr',-1.)
        call addpar('idtlgr',0.)

c |idtl|	The time step used for the release of particles. If this
c		is 0 particles are released only once at the beginning
c		of the lagrangian simulation. No particles are released 
c		for a value of less than 0. (Default 0)

        call addpar('idtl',0.)

c |itranf, itrend|	Initial and final time for the release of particles.
c			If not specified the particles are released over the
c			whole lagrangian simulation period.

        call addpar('itranf',-1.)
        call addpar('itrend',-1.)

c |ipvert|	Set the vertical distribution of particles:
c		\begin{description}
c		\item[0] releases one particles only in surface layer
c		\item[$>$0] release n particles regularly
c		\item[$<$0] release n particles randomly
c		\end{description}

        call addpar('ipvert',0.)

c |linbot|	Set the bottom layer for vertical releases (Default -1, bottom layer)

        call addpar('linbot',-1.)

c |lintop|	Set the top layer for vertical releases (Default 1, surface layer)
        call addpar('lintop',1.)

c |stkpar|	Calibration parameter for parameterizing the stokes drift 
c		induced by waves (and wind). Only affect particle of the sea
c		surface (layer = 1). The wind file is needed even in offline 
c		mode (Default 0). 

        call addpar('stkpar',0.)

c |dripar|	Parameter to account for drifter inertia by multiplying
c		the advective transports. Usually it assumes values between 
c		0.9 and 1.2 (Default 1). 

        call addpar('dripar',1.)

c |lbeach|	Parameter to account for particles beaching on the shore. It 
c		assumes values between 0 (no beaching) and 1 (Default 0).

        call addpar('lbeach',0.)

c |lgrlin|	File name that contains closed lines of the area where
c		the particles have to be released. If not given, the particles
c		are released over the whole domain.

        call addfnm('lgrlin',' ')

c |lgrtrj|	File name that contains a drifter trajectory with time 
c		(yyyy-mm-dd::HH:MM:SS) and position (x and y) of release 
c		of |nbdy| particles.

        call addfnm('lgrtrj',' ')

c |lgrini|	File name that contains initial particle distribution.
c		It has the same format as the lagrangian output.

        call addfnm('lgrini',' ')

c |itlgin|	Time to use for the initialization of the particle 
c		distribution from file (|lgrini|). 

        call addpar('itlgin',-1.)

cc To compute the transit time of the particles to leave
cc a specified area. |artype| is the flag to detect the
cc element outside the defined area. 
cc Default = -1, i.e., no transit times are computed

        call addpar('artype',-1.)

cc still to be commented

        call addpar('lcust',0.)
        call addpar('ldecay',0.)
        call addpar('ioil',0.)
        call addpar('ilarv',0.)
        call addpar('ised',0.)


c DOCS	END

	end

c************************************************************************

        subroutine nlsinh_wrt

        implicit none

c $wrt section

c DOCS  START   P_wrt
c
c Section |$wrt| contains parameters for computing water renewal time.
c During runtime if writes a .jas file with time series of total tracer
c concentration in the basin and WRT computed according to different methods.
c Nodal values of computed WRT are written in the .wrt file.
c Frequency distributions of WRTs are written in the .frq file.
c
c Please find all details here below.

        call sctpar('wrt')             !sets default section
        call sctfnm('wrt')

c |idtwrt|	Time step to reset concentration to c0. Use 0 if no reset
c		is desired. Use -1 if no renewal time computation is desired
c		(Default -1).

        call addpar('idtwrt',-1.)

c |itmin|	Time from when to compute renewal time (-1 for start of sim)
c		(Default -1)

        call addpar('itmin',-1.)

c |itmax|	Time up to when to compute renewal time (-1 for end of sim)
c		(Default -1).

        call addpar('itmax',-1.)

c |c0|		Initial concentration of tracer (Default 1).

        call addpar('c0',1.)

c |iaout|	Area code of elements out of lagoon (used for init and retflow).
c		Use -1 to if no outside areas exist. (Default -1).

        call addpar('iaout',-1.)

c |percmin|	Percentage to reach after which the computation is stopped.
c		Use 0 if no premature end is desired (Default 0).

        call addpar('percmin',0.)

c |iret|	Equal to 1 if return flow is used. If equal to 0 the 
c		concentrations outside are explicitly set to 0 (Default 1).

        call addpar('iret',1.)

c |istir|	If equal to 1 simulates completely stirred tank 
c		(replaces at every time step conz with average conz)
c		(Default 0).

        call addpar('istir',0.)

c |iadj|	Adjust renewal time for tail of distribution (Default 1).

        call addpar('iadj',1.)

c |ilog|	Use logarithmic regression to compute renewal time (Default 0).

        call addpar('ilog',0.)

c |ctop|	Maximum to be used for frequency curve (Default 0).

        call addpar('ctop',0.)

c |ccut|	Cut renewal time at this level (for res time computation)
c		(Default 0).

        call addpar('ccut',0.)

c |wrtrst|	If reset times are not regularly distributed (e.g., 1 month)
c		it is possible to give the exact times when a reset should
c		take place. |wrtrst| is a file name where these reset times
c		are specified, one for each line. For every line two integers
c		indicating date and time for the reset must be specified.
c		If only one value is given, time is taken to be 0. The format
c		of date is "YYYYMMDD" and for time "hhmmss". If the file
c		wrtrst is given |idtwrt| should be 0.

        call addfnm('wrtrst',' ')

c DOCS  END

        end

c************************************************************************

	subroutine nlsinh_bfmsc

c Parameters for BFM ECOLOGICAL module        !BFMSC

	implicit none

        call sctpar('bfmsc')             !sets default section
        call sctfnm('bfmsc')

c BFM ECOLOGICAL MODEL CALL 

	call addpar('ibfm',0.)
	call addpar('ibtanf',0.)
	call addpar('ibtend',0.)
	call addpar('itmbfm',-1.)
	call addpar('idtbfm',0.)
	call addpar('bligth',1.) !light flag=1 max/min light W/m**2 in nml file

	end

c************************************************************************

	subroutine nlsinh_proj

	implicit none

c $proj section

c DOCS  START   P_proj
c
c Section |$proj| handles the projection from cartesian to
c geographical coordinate system. If |iproj| $>$0 the projected geographical 
c coordinates can be used for computing spatially variable Coriolis parameter 
c and tidal potential even if the basin is in cartesian coordinate system 
c (|isphe| = 0) .
c
c Please find all details here below.

        call sctpar('proj')             !sets default section
        call sctfnm('proj')

c |iproj|	Switch that indicates the type of projection
c		(default 0):
c		\begin{description}
c		\item[0] do nothing
c		\item[1] Gauss-Boaga (GB)
c		\item[2] Universal Transverse Mercator (UTM)
c		\item[3] Equidistant cylindrical (CPP)
c		\item[4] UTM non standard
c		\end{description}

	call addpar('iproj',0.)

c |c\_fuse|	Fuse for GB (1 or 2, default 0)

	call addpar('c_fuse',0.)

c |c\_zone|	Zone for UTM (1-60, default 0)

	call addpar('c_zone',0.)

c |c\_lamb|	Central meridian for non-std UTM (default 0)

	call addpar('c_lamb',0.)

c |c\_x0|	x0 for GB and UTM (default 0)

	call addpar('c_x0',0.)

c |c\_y0|	y0 for GB and UTM (default 0)

	call addpar('c_y0',0.)

c |c\_skal|	Scale factor for non-std UTM (default 0.9996)

	call addpar('c_skal',0.9996)

c |c\_phi|	Central parallel for CPP (default 0.9996)

	call addpar('c_phi',0.)

c |c\_lon0|	Longitude origin for CPP (default 0)
	call addpar('c_lon0',0.)

c |c\_lat0|	Latitude origin for CPP (default 0)
	call addpar('c_lat0',0.)

c DOCS  END

	end

c************************************************************************

	subroutine nlsinh_undoc

c not documented parameters 

	implicit none

	call sctpar('para')		!sets default section

cc undocumented parameters

	call addpar('iclose',0.)
	call addpar('itsmed',0.)	!averages for T/S

cc internally used parameters

	call addpar('flag',0.)
	call addpar('zconst',-999.)	!set to flag
	call addpar('volmin',1.)	!minimum volume to remain in el.

cc debug

	call addpar('vreps',5.e-4)
	call addpar('vrerr',1.e-1)
	call addpar('levdbg',0.)	!debug level (0 = no, 9 = max)

cc distance for advective terms

	call addpar('nadist',0.)

cc new for scaling time step

	call addpar('itunit',1.)	!this is not supported anymore

cc experimental stuff

        !call addpar('nbsig',0.)         !sigma layers to read in for OBC

        call addpar('sedim',0.)         !sedimentation for theseus
        call addfnm('hsedim',' ')         !sedimentation hev file for theseus

        call addpar('nomp',0.)          !number of threads to use

	end

c************************************************************************

	subroutine nlsinh_georg

c parameters used for private projects

	implicit none

	call sctpar('para')		!sets default section

	call addpar('hlido',999.)	!maximum depth at lido
	call addpar('hmala',999.)	!maximum depth at lido
	call addpar('hchio',999.)	!maximum depth at lido

	call addpar('zrise',0.)		!sea level rise
	call addpar('zsalv',999.)	!save-guarding level
	call addpar('zfranc',0.)	!extra security for forecast

	end

c************************************************************************

	subroutine nlsinh_unused

c parameters not used anymore -> to be deleted

	implicit none

	call sctpar('para')		!sets default section

	call addpar('iprogr',0.)

	end

c ********************************************************************

c This subroutine defines the simulation wave parameter

        subroutine nlsinh_waves

        implicit none

c $waves section

c DOCS  START   P_waves
c
c Parameters in section |$waves| activate the wind wave module and define
c which kind of wind wave model has to be used. These parameters must
c be in section |waves|.

        call sctpar('waves')             !sets waves section
        call sctfnm('waves')

c |iwave|	Type of wind wave model and coupling procedure (default 0):
c		\begin{description}
c		\item[0] No wind wave model called 
c		\item[1] The parametric wind wave model is called 
c		(see file subwave.f)
c		\item[$>$1] The spectral wind wave model WWMIII is called
c		\item[2] ... wind from SHYFEM, radiation stress formulation
c		\item[3] ... wind from SHYFEM, vortex force formulation 
c		\item[4] ... wind from WWMIII, radiation stress formulation
c		\item[5] ... wind from WWMIII, vortex force formulation 
c		\item[11] The spectral wind wave model WaveWatch WW3 is called
c		\end{description}
c		When the vortex force formulation is chosen the wave-supported
c		surface stress is subtracted from the wind stress, in order to
c		avoid double counting of the wind forcing in the flow model.
c		Moreover, the use of the wave-depended wind drag coefficient 
c		could be adopted setting |itdrag| = 3.

        call addpar('iwave',0.)

c
c |dtwave|	Time step for coupling with WWMIII. Needed only for
c		|iwave| $>$ 1 (default 0).

        call addpar('dtwave',0.)

c |idtwav|, |itmwav|	Time step and start time for writing to file wav,
c			the files containing output wave variables (significant 
c			wave height, wave period, mean wave direction).

        call addpar('idtwav',0.)
        call addpar('itmwav',-1.)

c
c DOCS  END
c

        end

c************************************************************************

	subroutine nlsinh_nonhydro

c parameters for non hydrostatic model (experimental)

	implicit none

        call sctpar('nonhyd')           !sets default section

	call addpar('inohyd',0.)	! 0=hydrostatic  1=non-hydrostatic

        call addpar('aqpar',0.5)	!implicit parameter

	call addpar('islope',0.)	!type of grid for poisson equation

        call addpar('ivwadv',0.)        !vert advection of vert momentum
        call addpar('inhflx',0.)        !flux upwind for horiz advect of w
	call addpar('inhadj',0.)        !choice for correction of U,V,eta
        call addpar('inhwrt',0.)        !output every inhwrt time steps
        call addpar('inhbnd',0.)        !exclude NH dynamics for boundaries
        call addpar('iwvel',0.)         !write vertical velocity
        call addpar('iqpnv',0.)         !write NH pressure !DWNH
        call addpar('nqdist',0.)        !distance for NH pressure terms

	end

c************************************************************************

	subroutine nlsinh_connect

c parameters for connectivity (experimental)

	implicit none

        call sctpar('connec')          !sets default section
        call sctfnm('connec')

	call addpar('icnn',0.)	        !number of stations for connectivity

        call addpar('radcnn',0.)	!radius of release area
        call addpar('ppscnn',0.)	!particles per second to be released
        call addpar('pldcnn',0.)	!pelagic larval duration
        call addpar('idtcnn',0.)	!output every idtcnn

        call addfnm('statcnn',' ')      !connectivity file with stations

	end

c************************************************************************
c************************************************************************
c************************************************************************
c************************************************************************
c************************************************************************

	subroutine nlsina

c initializes the parameter file for the post processing routines

	implicit none

	call nlsina_general
	call nlsina_para
	call nlsina_color
	call nlsina_arrow
	call nlsina_legend
	call nlsina_legvar
	call nlsina_sect

	end

c************************************************************************

	subroutine nlsina_general

c $para section with general variables

	implicit none

c $para section for plotting

c DOCS	START	S_para_general
c
c The next parameters decide what to plot. They can be set directly in
c the STR file here. However, the preferred way to set them is through
c the back end |plots| that automatically sets these parameters.

	call sctpar('para')		!sets default section
	call sctfnm('para')		!sets default section

c The parameter |iwhat| is compulsory and defines the variable to be plotted.
c Please note again that this parameter will be set if invoking the plot
c program through |plots|.
c
c |iwhat|		Flag that determines what to plot. If 0 then
c			the program asks interactively for it. (Default 0)
c			\begin{description}
c			\item[1] Plot basin (grid and isolines of depth)
c			\item[2] Plot velocities
c			\item[3] Plot transports
c			\item[4] Plot water levels
c			\item[5] Plot concentration
c			\item[6] Plot temperature
c			\item[7] Plot salinity
c			\item[8] Plot rms-velocity
c			\item[9] (not used)
c			\item[10] Plot generic scalar (see |ivar|)
c			\item[11] Plot wind vectors
c			\item[12] Plot lagrangian particles
c			\item[13] Plot wave data
c			\end{description}

	call addpar('iwhat',0.)		!what to plot

c The next parameters define other aspects of the plot. In particular their
c meaning is
c
c |iauto|		Normally the simulation name and basin are
c			asked interactively when running the plotting
c			program. However, if |iauto| is different from 0
c			then the file |.memory| is read and the information
c			contained in it is used. The file memory can be
c			set through the program |memory| and it is
c			changed when other parameters are inputted
c			interactively.
c |level|		For 3D applications it indicates the vertical level
c			for which the plot is desired. 1 indicates the 
c			surface layer, 2 the second layer from the surface etc.
c			0 gives integrated quantities, and a value of -1
c			indicates the bottom layer (which refers to different
c			layers for every element). (Default 0)
c |ivar|		Variable to be plotted for scalar quantities. In the
c			file NOS more then one scalars can be contained. To
c			choose the desired scalar, the id of the scalar has 
c			to be given in |ivar|. Note that if only one variable
c			type is contained in the file, then it is not
c			necessary to set |ivar|.

	call addpar('iauto',0.)		!silent mode
	call addpar('level',0.) 	!level (-1 -> bottom   0 -> integr.)
	call addpar('ivar',0.)		!what variable to plot

cc still to document

	call addfnm('varnam',' ')	!variable name to be plotted
	call addpar('isect',0.)		!vertical section

c The next variables define the time of the plots. Even if the names of two
c of the variables are identical to variables used in the finite element
c model, the meaning is different. In the simulation model, |itanf, itend|
c define the initial and final time of simulation, whereas here these
c variables define the initial and final time of the plot of results.
c if none of the time variables are set, then all time records are plotted.

c |itanf|		Initial time for plotting. (Default is 
c			first data record)
c |itend|		Final time for plotting. (Default is 
c			last data record)
c |nout|		Frequency for plotting. A value of 0 or 1
c			plots every record, 2 plots every other record etc.
c			A negative value works as a positive one,
c			but starts plotting from the first record. So
c			-2 plots records 1,3,5,etc.
c			(Default 1) 

	call addpar('itanf',-1.)	!time start
	call addpar('itend',-1.)	!time end
	call addpar('nout',1.)		!time frequency

	call addpar('atanf',-1.)	!time start (absolute)
	call addpar('atend',-1.)	!time end (absolute)

c DOCS	END

	end

c************************************************************************

	subroutine nlsina_para

c $para section

	implicit none

c DOCS	START	S_para_a
c
c The parameters in section |$para| set generic values for the plot.

	call sctpar('para')		!sets default section
	call sctfnm('para')		!sets default section

c Some of the parameters set coordinates in the plot. For example, the
c values |x0, y0| and |x1, y1| indicate the actual plotting area, which can
c be bigger or smaller than the extension of the numerical grid.
c
c Normally, values have to be in meters (the same as the coordinates in the
c numerical grid). However, also relative coordinates can be used. If all
c values given are in the range between -1 and +2, then these values
c are interpreted as relative coordinates. Therefore, $x$ coordinates of
c 0 indicate the left border, and 1 the right border. The upper left quarter
c of the domain can be chosen with (|x0, y0|) = (0,0.5) and
c (|x1, y1|) = (0.5,1).

c |x0, y0|		Lower left corner of the plotting area.
c			(Default is whole area)
c |x1, y1|		Upper right corner of the plotting area.
c			(Default is whole area)

	call addpar('x0',0.)		!dimension of plot
	call addpar('y0',0.)		!dimension of plot
	call addpar('x1',0.)		!dimension of plot
	call addpar('y1',0.)		!dimension of plot

c The next values give the position, where the legend (scale bar and
c true north) is plotted. This legend will only be plotted if
c the coordinates are not geographical (lat/lon) but cartesian.

c |x0leg, y0leg|	Lower left corner of the area
c			where the legend is plotted.
c |x1leg, y1leg|	Upper right corner of the area.
c			where the legend (north and scale) is plotted.

	call addpar('x0leg',0.)		!dimension of legend
	call addpar('y0leg',0.)		!dimension of legend
	call addpar('x1leg',0.)		!dimension of legend
	call addpar('y1leg',0.)		!dimension of legend

c |lblank|		The legend is plotted over a white rectangle.
c			Sometimes this blanking is not desirable.
c			If you do not want to have a white box below the
c			legend set |lblank| to 0. (Default 1)

	call addpar('lblank',1.)

c |cislnd|		It is possible to plot all islands in gray color.
c			Setting |cislnd| to a value between 0 (black) and 
c			1 (white) will achieve this. A negative value
c			will not fill islands with gray color.
c			(Default -1)

	call addpar('cislnd',-1.)	!plot also outer island: 2 <= c <= 3

c |dgray|		It is possible to plot all dry areas in gray color.
c			Setting |dgray| to a value between 0 (black) and 
c			1 (white) will achieve this. A negative value
c			will not fill dry areas with gray color.
c			(Default -1)

c |hgray|		Whereas |dgray| is normally only coloring
c			elements that are dry, you can also color elements
c			shallower than a given depth |hgray|. E.g., a value
c			for |hgray| of -0.5 will plot in gray all
c			elements with depth lower than -0.5 m (salt
c			marshes). (Default -10000)

	call addpar('dgray',-1.)
	call addpar('hgray',-10000.)	!gray all elems with h < hgray

c |dxygrd|		Grid size if the results are interpolated on
c			a regular grid. A value of 0 does
c			not use a regular grid but the original
c			finite element grid for plotting. (Default 0)
c |typls|		Typical length scale to be used when scaling
c			velocity or transport arrows. If |dxygrd| is
c			given this length is used and |typls| is not used.
c			If not given it is computed from the basin
c			parameters. (Default 0)
c |typlsf|		Additional factor to be used with typls to
c			determine the length of the maximum or
c			reference vector. This is the easiest way
c			to scale the velocity arrows 
c			with an overall factor. (Default 1)
c |velref|		Reference value to be used when scaling arrows.
c			If given, a vector with this value will have a length
c			of |typls|*|typlsf| on the map, or, in case
c			|dxygrd| is given, |dxygrd|*|typlsf|. If not set
c			the maximum value of the velocity/transport
c			will be used as |velref|. (Default 0)
c |velmin|		Minimum value for which an arrow will be plotted.
c			With this value you can eliminate small arrows
c			in low dynamic areas. (Default 0)
c |velmax|		Maximum value for which an arrow will be plotted.
c			With this value you can eliminate arrows that are
c			too big. This is useful if you would like to study an
c			area with low current speed but adjacent area have
c			high current speeds that would overplot the area.
c			(Default -1, no limitation)

	call addpar('dxygrd',0.)	!grid size for regular grid
	call addpar('typls',0.)		!typical length scale for arrow plot
	call addpar('typlsf',1.)	!factor for typical length scale
	call addpar('velref',0.)	!reference velocity for length scale
	call addpar('velmin',0.)	!minimum velocity to be plotted
	call addpar('velmax',-1.)	!maximum velocity to be plotted

c |isphe|		If 0 a cartesian coordinate system is used,
c			If 1 the coordinates are in the spherical 
c			system (lat/lon). Among other, this
c			indicates that the $x$-coordinates will be multiplied
c			by a factor that accounts for the visual deformation
c			using lat/lon coordinates.
c			The default is -1, which means that the 
c			type of coordinate system will 
c			be determined automatically. (Default -1)
c |reggrd|		If different from 0 it plots a regular grid over
c			the plot for geographical reference. The value of
c			|reggrd| gives the spacing of the regular grid lines.
c			The units must be according to the units used for
c			the coordinates. With value of -1 the regular grid is
c			determined automatically. (Default -1)
c |regdst|		This value gives the number of intervals
c			that are used to sub-divide the grid given by
c			|reggrd| with a black and white scale around
c			the plot. If 0 it tries to determine automatically
c			the sub-intervals (2 or 4). A value of -1 does
c			not plot the subgrid scale. (Default 0)
c |reggry|		If plotting the regular overlay grid this gives
c			the gray value used for the grid. 0 is black, and
c			1 is white. A value of 1 does not plot the
c			overlay grid, but still writes the labels. 
c			(Default 1)

	call addpar('isphe',-1.)	!spherical coordinate system
	call addpar('reggrd',-1.)	!regular grid spacing
	call addpar('regdst',0.)	!regular micro grid spacing
	call addpar('reggry',1.)	!gray value

c |bndlin|		Name of file that gives the boundary line
c			that is not part of the finite element domain.
c			The file must be in GRD format. An older BND
c			format is also accepted, but deprecated.
c			(Default is no file)

	call addfnm('bndlin'," ")	!name of boundary line file

c |ioverl|		Create overlay of velocity vectors on scalar value.
c			With the value of 0 no overlay is created, 1
c			creates an overlay with the velocity speed.
c			The value of 2 overlays vertical velocities
c			3 water levels and 4 overlays bathymetry.(Default 0)
c |inorm|		Normally the horizontal velocities are plotted
c			in scale. The value of |inorm| can change this
c			behavior. A value of 1 normalizes velocity vectors
c			(all vectors are the same length), whereas 2
c			scales from a given minimum velocity |velmin|.
c			Finally, the value of 3 uses a logarithmic scale.
c			(Default 0)

	call addpar('ioverl',0.)	!overlay in color
	call addpar('inorm',0.)		!vertical velocity as overlay

c The next parameters give the choice to selectively avoid to plot areas
c of the basin and to apply different gray tones for the boundary and
c net lines when plotting the basin.
c Please remember that when working with gray tones the value should
c be between 0 (black) and 1 (white).
c
c |ianopl|	Area code for which no plot has to be produced. Normally 
c		the whole basin is plotted, but with this parameter some
c		areas can be excluded. (Default -1)
c |bgray|	Gray value used for the finite element grid when plotting
c		the bathymetry. (Default 0.8)
c |bbgray|	Gray value used for the boundary of the finite element grid.
c		(Default 0)
c |bsgray|	Gray value used to plot the finite element grid over
c		a scalar or velocity plot. This is basically useful
c		for debugging reasons. The default is to not plot
c		the grid (Default -1.0)

        call addpar('ianopl',-1.)      !do not plot these areas
	call addpar('bgray',0.8)       !gray value for bathymetry
	call addpar('bbgray',0.0)      !gray value for boundary
	call addpar('bsgray',-1.0)     !gray value for plotting maps

c The next two parameters handle the plotting of the lagrangian particles.
c
c |lgrtrj|	If equal 1 plot trajectories
c		instead of particle position. (Default 0)

	call addpar('lgrtrj',0.)	!lagrangian trajectories or position

c |lgmean|	Plot mean positions/trajectories.
c		With the value of 0 no mean pos/traj are created, 1
c		plot mean pos/traj together with single values, 
c		2 plot only mean pos/trajs, 3 as 2 but the first
c		trajectory is plot in tichk line. (Default 0)

	call addpar('lgmean',0.)	!lagrangian mean pos/trajs

c The next parameters handle the plotting of the basin.

c |ibox|	If set to 1 plots box information if the area code is
c		used to indicate the box number. This parameter
c		is only useful for the plotting of the box model. (Default 0)
c |inumber|	If set to 1 plots node and element numbers on top
c		of the grid. This is only useful in debug mode. (Default 0)

        call addpar('ibox',0.)
        call addpar('inumber',0.)

c DOCS	END

cc not documented

	call addpar('ifreg',0.)		!plot regular grid from fem file
	call addpar('iexreg',1.)	!plot half valid boxes in reg grid

	call addfnm('obspnt'," ")	!name of file with obs points
	call addfnm('metpnt'," ")	!name of file with meteo points
	call addfnm('spcvel'," ")	!name of file for velocity points

cc only for compatibility ... are not used anymore

	call addpar('traref',0.)	!reference transport for length scale
	call addpar('tramin',0.)	!minimum transport to be plotted

cc internally needed

	call addpar('dirn',0.)
	call addpar('href',0.)
	call addpar('hzmin',0.01)
	!call addpar('hzoff',0.05)
	!call addpar('hlvmin',0.25)

	end

c************************************************************************

	subroutine nlsina_color

c $color section

	use para

	implicit none

c DOCS	START	S_color
c
c Section |$color| deals with the definition of the colors
c to be used in the plot. A color bar is plotted too.

	call sctpar('color')		!sets default section
	call sctfnm('color')		!sets default section

c |icolor|	Flag that determines the type of color table
c		to be used. 0 stands for gray scale, 1 for
c		HSB color table. Other possible values are
c		2 (from white to blue), 3 (from white to red),
c		4 (from blue over white to red) and 
c		5 (from blue over black to red).
c		Values 6 and 7 indicate non-linear HSB color tables.
c		(Default 0)

	call addpar('icolor',0.)	!use color (1.) or not (0.)

c |colfil|	A color table can also be read from file. An example
c		of the format can be found in directory |femplot/color|
c		in the file |colormap.dat|. The variable |colfil|
c		indicates the file where the color table is being
c		read from. The default is not to read a color table file.
c |coltab|	If a color table file has been read then the variable
c		|coltab| indicates the name of the color table that
c		is going to be used. The default is to not use any
c		of the color tables if no name is specified.

        call addfnm('colfil',' ')
        call addfnm('coltab',' ')

c |isoval|	Array that defines the values for the isolines
c		and colors that are to be plotted. Values given must
c		be in the unit of the variable that will be plotted,
c		i.e., meters for water levels etc. 
c |color|	Array that gives the color indices for the
c		plotting color to be used. Ranges are from
c		0 to 1. The type of the color depends on the 
c		variable |icolor|. For the gray scale table
c		0 represents black and 1 white. Values in between
c		correspond to tones of gray. For the HSB color table
c		going from 0 to 1 gives the color of the rainbow.
c		There must be one more value in |color| than in |isoval|.
c		The first color in |color| refers to values less
c		than |isoval(1)|, the second color in |color| to
c		values between |isoval(1)| and |isoval(2)|. The last
c		color in |color| refers to values greater than the last
c		value in |isoval|.

	!call addpar('isoval',0.)
	!call addpar('color',0.)
	call para_add_array_value('isoval',0.)
	call para_add_array_value('color',0.)

c |x0col, y0col|	Lower left corner of the area where the
c			color bar is plotted.
c |x1col, y1col|	Upper right corner of the area where the
c			color bar is plotted.

	call addpar('x0col',0.)		!dimension of color bar
	call addpar('y0col',0.)		!dimension of color bar
	call addpar('x1col',0.)		!dimension of color bar
	call addpar('y1col',0.)		!dimension of color bar

c |cblank|		The color bar is plotted over a white rectangle.
c			Sometimes this blanking is not desirable.
c			If you do not want to have a white box below the
c			legend set |cblank| to 0. (Default 1)

	call addpar('cblank',1.)

c |faccol|	Factor for the values that are written to the 
c		color bar legend. This enables you, e.g., to give water level
c		results in mm (|faccol = 1000|). (Default 1)
c |ndccol|	Decimals after the decimal point for the values
c		written to the color bar legend. Use the value |-1|
c		to not write the decimal point. A value of 0 automatically
c		computes the number of decimals needed. (Default 0)
c |legcol|	Text for the description of the color bar. This text
c		is written above the color bar.

	call addpar('faccol',1.)	!factor for color bar
	call addpar('ndccol',-1.)	!decimals after point
	call addfnm('legcol'," ")	!legend for colorbar

c It is not necessary to give all values for isolines and colors above.
c A faster way is to give only the minimum and maximum values and fix
c the number of isovalues to be used.

c |niso|		Total number of isolines to use. (Default is |nisodf|)
c |nisodf|		Default number of isolines to use. (Default 5)
c |colmin, colmax|	Minimum and maximum color index used. Defaults are
c			0.1 and 0.9 respectively. The value of |colmax| can
c			be smaller than |colmin| which inverts the color
c			index used.
c |valmin, valmax|	Minimum and maximum value for isovalues to be used.
c			There is no default.
c |rfiso|		Defines function to be used to compute intermediate
c			values between |valmin| and |valmax|. If 0 or 1
c			the values are linearily interpolated. Else they
c			are computed by $y=x^n$ where $n$ is |rfiso|
c			and $x=\frac{v-v_{min}}{v_{max}-v{min}}$. Values
c			for |rfiso| greater than 0 capture higher detail in
c			the lower values, whereas values less than 1 do
c			the opposite.
c			(Default 0)
c |ipllog|		Indicates the usage of a logarithmic color scale.
c			The possible values are 0-3. The value of 0
c			indicates not to use a logarithmic scale.
c			If 1, the values of
c			the scale are 1,10,100,etc., if 2 the values
c			1,2,10,20,100,etc. are used, and for 3 the values
c			are 1,2,5,10,20,50,100,etc. (Default 0)
c |dval|		Difference of values between isolines. If this
c			value is greater then 0 the values for isolines 
c			and the total number of isolines are computed 
c			automatically using also |valmin| and |valmax|. 
c			(Default 0)

        call addpar('niso',0.)         !total number of isovalues
        call addpar('nisodf',5.)       !default number of isovalues to use
        call addpar('colmin',0.1)      !min color [0..1]
        call addpar('colmax',0.9)      !max color [0..1]
        call addpar('valmin',0.)       !min isovalue
        call addpar('valmax',0.)       !max isovalue
	call addpar('rfiso',0.)	       !function for intermediate values
	call addpar('ipllog',0.)       !logarithmic scale
	call addpar('dval',0.)	       !increment for autom. color sep.

c Since there is a great choice of combinations between the parameters,
c we give here the following rules how the values for colors and isolines
c are determined.
c
c If colors are given in array |color|, they are used, else |colmin| and
c |colmax| or their respective defaults are used to determine the color bar.
c If |isoval| is given it is used, else |valmin| and |valmax| are used.
c If |valmin| and |valmax| are not given they are computed every time
c for each plot and the minimum and maximum value in the basin are used.
c In any case, if |isoval| is specified the total number of isovalues
c is known and |niso| is ignored. However, if |isoval| is not given
c then first |dval| is used to decide how many isovalues to plot, and
c if |dval| is 0 then the |niso| and finally |nisodf| is used.
c
c Other parameters that can be changed are the following.

c |nisomx|	Maximum for |niso| allowed. This is especially useful
c		when the value for |niso| is determined automatically.
c		It avoids you to plot 1000 isolines due to wrong settings
c		of |dval|. However, if you want to use 50 isovalues
c		then just set |niso| and |nisomx| to 50. (Default 20)
c |nctick|	Number of values to be written in color bar. If |niso| is high
c		the labels on the color bar become unreadable. Therefore
c		you can use |nctick| to write only some of the values to
c		the color bar. For example, if |valmin| is 0 and |valmax| is
c		5 and you use many isolines, then setting |nctick| to 6 would
c		give you labels at values 0,1,2,3,4,5. If |nctick| is 0
c		then all labels are written. (Default 0)
c |isolin|	Normally the isolines are not drawn on the plot, just
c		the colors are used to show the value in the different
c		parts of the plot. A value different from 0 plots also 
c		the isolines. In this case |isolin| gives the number of
c		isolines to be plotted. A good choice is to make this
c		equal to |nctick|, so that the isolines correspond to the 
c		values	written on the colorbar. For compatibility, a value of
c		1 plots all isolines. (Default 0)
c |isoinp|	Normally inside elements the values are interpolated.
c		Sometimes it is useful to just plot the value of the
c		node without interpolation inside the element. This can
c		be accomplished by setting |isoinp=0|. Setting instead
c		|isoinp| to a value of 2 plots a constant value in
c		the element. (Default 1)

        call addpar('nisomx',20.)      !maximum number of isovalues allowed
        call addpar('nctick',0.)       !default number of ticks to use
        call addpar('isolin',0.)       !plot isolines with color ?
        call addpar('isoinp',1.)       !interpolate inside elements

c DOCS	END

	end

c************************************************************************

	subroutine nlsina_arrow

c $arrow section

	implicit none

c DOCS	START	S_arrow
c
c Parameters in section |$arrow| deal with the reference arrow that is plotted
c in a legend. The arrow regards the plots where the velocity or
c the transport is plotted.

	call sctpar('arrow')		!sets default section
	call sctfnm('arrow')		!sets default section

c |x0arr, y0arr|	Lower left corner of the area where the
c			reference arrow is plotted.
c |x1arr, y1arr|	Upper right corner of the area where the
c			reference arrow is plotted.

	call addpar('x0arr',0.)		!dimension of color bar
	call addpar('y0arr',0.)		!dimension of color bar
	call addpar('x1arr',0.)		!dimension of color bar
	call addpar('y1arr',0.)		!dimension of color bar

c |ablank|		The arrow legend is plotted over a white rectangle.
c			Sometimes this blanking is not desirable.
c			If you do not want to have a white box below the
c			legend set |ablank| to 0. (Default 1)

	call addpar('ablank',1.)

c |facvel|		Factor for the value that is written to the 
c			arrow legend for the velocity.
c			This enables you, e.g., to give 
c			velocities in mm/s (|facvel = 1000|). (Default 1)
c |ndcvel|		Decimals after the decimal point for the values
c			written to the arrow legend. 
c			Use the value |-1|
c			to not write the decimal point. (Default 2)
c |legvel|		Text for the description of the arrow legend.
c			This text is written above the arrow.
c |arrvel|		Length of arrow in legend (in velocity
c			units). If not given the arrow length will be computed
c			automatically. (Default 0)
c |sclvel|		Additional factor to be used for the arrow
c			in the legend. When the arrow length will be
c			computed automatically, this parameter gives
c			the possibility to change the length of the
c			reference vector. This is an easy way
c			to scale the velocity arrow
c			with an overall factor. Not used if
c			|arrvel| is given. (Default 1)

	call addpar('facvel',1.)	!factor for velocity
	call addpar('ndcvel',2.)	!decimals after point (velocity)
	call addfnm('legvel'," ")	!legend for velocity
	call addpar('arrvel',0.)	!length of arrow
	call addpar('sclvel',1.)	!factor for arrow

c DOCS	END

	end

c************************************************************************

	subroutine nlsina_legend

c $legend section

	implicit none

c DOCS	START	S_legend

c In section |$legend| annotations in the plots can be given. The
c section consists of a series of lines that must contain the 
c following information:
c
c The first value is a keyword that specifies what has to be plotted.
c Possible values are |text|, |line|, |vect|, |rect|, |circ|  and also
c |wid| and |col|. These correspond to different types of information
c that is inserted into the plot such as text, line, vector, rectangle
c or circle (filled or just outline). Moreover, the color and
c line width of the pen can be controlled by with |wid| and |col|.
c
c In case of |text| the starting position (lower left corner) is given,
c then the point size of the font and the text that is inserted. |line|
c needs the starting and end point of the line. The same with |vect|,
c but in this case also the relative tip size must be given as a final
c parameter. |rect| needs the coordinates of the lower left corner and 
c upper right corner of the rectangle. It also needs the color used for
c the filling of the rectangle (0-1) or the flag -1 which only draws the
c outline of the rectangle without filling it. |circ| needs the center
c point and the radius and a fill color (see rectangle). Finally |wid| needs 
c the relative width of the line and |col| the stroke color used when plotting
c lines.
c
c A small example of an annotation that explains the above parameters
c would be:
c
c \begin{verbatim}
c $legend
c text  30500 11800     15  'Chioggia'   #text, 15pt
c line  30500 11800 35000 15000          #line
c vect  30500 11800 35000 15000 0.1      #arrow, tip size 0.1
c rect  30500 11800 35000 15000 0.1      #rectangle, fill color 0.1
c rect  30500 11800 35000 15000 -1       #rectangle (outline, no fill)
c circ  30500 11800 5000 -1              #circle (outline, no fill)
c wid   5                                #set line width to 5
c col   0.5                              #set color to 0.5
c $end
c \end{verbatim}
c
c There is also an old way to specify the legend that does not use
c keywords. However, this way is deprecated and unsupported and is therefore
c not described anymore in this manual.

c DOCS	END

	end

c************************************************************************

	subroutine nlsina_legvar

c $legvar section (variable legend)

	implicit none

c DOCS	START	S_legvar
c
c In section |$legvar| variable fields like the date and wind vectors
c may be inserted into the plot.

	call sctpar('legvar')		!sets default section
	call sctfnm('legvar')		!sets default section

c A time and date can be assigned to the simulation results. These values
c refer to the time 0 of the FEM model. The format for the date is
c YYYYMMDD and for the time HHMMSS. 
cc Please note that the date should not be
cc given as YYYYMMDD because due to precision problems this will not work. 
c You can also give a time zone if your time is not referring to 
c GMT but to another time zone such as MET. Please note that you have to give
c this information only if the simulation does not contain it already.
c Normally, this information is already assigned during the simulation runs.

c |date|		The real date corresponding to time 0. (Default 0)
c |time|		The real time corresponding to time 0. (Default 0)
c |tz|			The time zone you are in. This is 0 for GMT, 1 for MET
c			and 2 for MEST (MET summer time). (Default 0)
c |tzshow|		The time zone you want to show your results. If your
c			time zone is GMT (0) and you want to show the results
c			referred to MET (+1) set this to +1. Please note that
c			you have to set this variable only if you want to
c			show results in a different time zone than the one
c			given in |tz|. (Default 0)

	call addpar('date',-1.)
	call addpar('time',0.)
	call addpar('tz',0.)
	call addpar('tzshow',0.)

c The information of date and time may be written to the plot. This
c is done with the following parameters.

c |xdate, ydate|	Starting point for the date text (lower left corner).
c |sdate|		Point size of the text. (Default 18)
c |idate|		Output mode. If 0 no date is written to the
c			plot, else the date and time is written. (Default 0)

	call addpar('xdate',0.)
	call addpar('ydate',0.)
	call addpar('sdate',18.)         !size
	call addpar('idate',0.)         !mode

c Wind data can be used to insert a wind vector into the figure.
c This is useful because in the case of variable wind 
c the direction and speed of the wind that was blowing
c in the moment of the plot is shown.
c
c Since only one wind vector can be plotted, the wind data must consist
c of one value for each time. The same ASCII file that is used
c in the STR file can be used.

c |xwind, ywind|	Starting point where the wind arrow is plotted.
c |iwtype|		Type of wind data. The same as the one in the
c			STR file. If this parameter is 0 then no
c			wind vector is plotted. (Default 0)
c |lwwind|		Line width of the wind vector. (Default 0.1)
c |scwind|		Scaling parameter of the wind vector. This depends
c			on the size of your plot. If your wind is 10 m/s
c			and you want the vector to stretch over a distance
c			of 5 km on the plot then you have to choose
c			the value of 500 (10*500=5000) for |scwind|.
c			(Default 1)
c |wfile|		Name of the file containing the wind data. This 
c			may be the same file than the one used in the
c			STR file to run the program.

	call addpar('xwind',0.)
	call addpar('ywind',0.)
	call addpar('iwtype',1.)         !mode
	call addpar('lwwind',0.1)        !line width
	call addpar('scwind',1.)         !scale
	call addfnm('wfile',' ')	 !name of wind file

c The wind vector is also given a text legend with the speed of the
c wind written out. The next parameters decide where and how this information
c is put.

c |xtwind, ytwind|	Starting point for the legend text (lower left corner).
c |stwind|		Point size of the text. (Default 18)
c |wtext|		Text used for the legend (Default 'Wind speed')
c |wunit|		Unit for the wind speed (Default 'm/s')

	call addpar('xtwind',0.)
	call addpar('ytwind',0.)
	call addpar('stwind',18.)         !size
	call addfnm('wtext','Wind speed') !legend for wind
	call addfnm('wunit','m/s')	  !unit for wind

c DOCS	END

	end

c*******************************************************************

	subroutine nlsina_sect

c $sect section (vertical section)

	implicit none

c DOCS	START	S_sect
c
c In this section a vertical section plot can be set up. This will
c allow to plot 3D variables not horizontally for each layer, but
c along a given section in the vertical

	call sctpar('sect')		!sets default section
	call sctfnm('sect')		!sets default section

c |vsect|		Name of the file containing the node list defining
c			the section. The nodes must be adjacent in the
c			numerical grid used. There should be one node
c			on each line. The program |make_line_nodes.sh|
c			can be used to produce such a node list from a
c			given line created with |grid|. If this file
c			is not given, no vertical section will be plotted.

	call addfnm('vsect','')		!name of line defining vertical section

c The next parameters decide about how the plot is scaled in the vertical
c and if a grid overlay is created.
c
c |ivert|		This parameter decides what is plotted on the
c			vertical axis. If 0, depth is plotted on the axis.
c			If 1, vertical layers are plotted. If 2, depth is
c			plotted, but it is scaled with a logarithm, giving
c			more importance to the layers close to the surface.
c			(Default 0)
c |ivgrid|		If 1 a grid is plotted as overlay over the
c			vertical section plot. This is mainly needed
c			for a better orientation where the nodes and 
c			depth values or layers are situated. (Default 0)
c |hvmax, lvmax|	Normally the whole vertical extension is plotted
c			for the section. Sometimes only the upper part
c			of the water column may be needed. With these
c			two parameters the maximum depth (|hvmax|) or
c			the maximum layer (|lvmax|) to be plotted can
c			be imposed. The values may be also higher than
c			the maximum depth/layer. In this case extra space
c			is added below the plotting area. Only one
c			of the two parameters can be set.
c |bsmt|		Set to 1 to have smooth bottom profile even in 
c			case of z-layer structure. (Default 0)

	call addpar('ivert',0.)		!0: depth  1: layers  2: log depth
	call addpar('ivgrid',0.)	!plot grid layout over plot
	call addpar('hvmax',0.)		!maximum depth to be plotted
	call addpar('lvmax',0.)		!maximum layer to be plotted
	call addpar('bsmt',0.)		!1: smooth bottom

c The vertical section plot also creates a color bar. This color bar is
c normally put on the right side of the plot. If there is space inside the
c plot it might be placed on top of the section plot. In this case the
c following parameters can be used. Please note that in this case it is
c the relative position (from 0 to 1) that has to be specified.
c
c |x0sect, y0sect|	Lower left corner of the area where the
c			color bar is plotted.
c |x1sect, y1sect|	Upper right corner of the area where the
c			color bar is plotted.

	call addpar('x0sect',0.)
	call addpar('y0sect',0.)
	call addpar('x1sect',0.)
	call addpar('y1sect',0.)

c The following parameters set titles for the plot, the axis, and for
c an extra description of the starting and end point of the plot.
c They have some (intelligent) default values.
c
c |vtitle|		Title of the plot.
c |xtitle|		Title for the x-axis.
c |ytitle|		Title for the y-axis.
c |ltitle|		Title for start point of the line. (No default)
c |rtitle|		Title for end point of the line. (No default)

	call addfnm('vtitle','Section Plot')	!title for plot
	call addfnm('xtitle','Distance [m]')	!title for x-axis
	call addfnm('ytitle','Depth [m]')	!title for y-axis
	call addfnm('ltitle','')		!title for left point
	call addfnm('rtitle','')		!title for right point

c When plotting velocities you can decide if using for the color the normal 
c velocity across the section or the tangent velocity.
c In any case, in order to visualize also the velocity
c tangent to the section arrwos are used. The next parameters deal with
c the scaling of these arrows.

c |vmode|	When plotting velocities as default the normal velocity 
c		across the section is used for the color plot (|vmode|=0).
c		If you want to use the tangential velocity, please set
c		|vmode|=1. (Default 0)
c |avscal|	This parameter defines the horizontal scale for the 
c		velocity vector. It defines the length scale 
c		in units of x-coordinates of the vertical section so
c		that a velocity of 1 m/s fits comfortably
c		into the plot. If 0 the scale is computed automatically.
c		Please note that in this case the velocities will be
c		scaled differently for every plot. Setting |avscal| 
c		therefore guarantees that the velocity arrows will
c		be comparable between plots. If |avscal| is negative,
c		instead of using x-coordinates units, the legth scale is
c		in [cm]. Therefore, a value of -1 will represent a
c		velocity of 1 m/s with 1 cm on the plot. This unit for
c		the length scale is sometimes easier to understand. 
c		(Default 0)
c |rvscal|	Extra factor that multiplies the scale factor. If your
c		automatic scale gives you vectors which are too small, you
c		can use |rvscal| to increase them. (Default 1)
c |rwscal|	Extra factor for the vertical scale. Normally the
c		vertical scale is computed automatically, If you dont
c		like the size of the vertical vectors you can controll
c		it with this parameter. A value of 2 will give you
c		a vertical vector twice as big a the default. (Default 1)
c |dxmin|	Sometimes there are two many arrows plotted horizontally.
c		The value of |dxmin| gives the minimum distance arrows have
c		to be apart to be plotted. A value of 0 plots all arrows.
c		(Default 0)
c |svtip|	The (relative) tip size of the arrow. It specifies how
c		big the arrow will be drawn. A value of 0 only draws the
c		arrow line without tip, and a negative value inhibits
c		drawing arrows at all. (Default 0.3)
c |rxscal,ryscal|	In case arrows are plotted, also a reference
c			vector is plotted. The size of this reference
c			vector s computed automatically, but can be 
c			controlled additionally by the parameters
c			|rxscal,ryscal|, which are in relative units
c			with respect to the reference box plotted.
c			(Default 0.6)

	call addpar('vmode',0.)
	call addpar('avscal',0.)
	call addpar('rvscal',1.)
	call addpar('rwscal',1.)
	call addpar('dxmin',0.)
	call addpar('svtip',0.3)
	call addpar('rxscal',0.6)
	call addpar('ryscal',0.6)

cc still to do: plot vector legend

c DOCS	END

	end

c*******************************************************************
c*******************************************************************
c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine fnminh

c initializes default names

	use para

	implicit none

	logical haspar

c $name section

	call sctfnm('name')		!sets default section

	call para_deprecate('basdir','')
	call para_deprecate('datdir','')
	call para_deprecate('tmpdir','')
	call para_deprecate('defdir','')

c DOCS	START	S_name
c
c In section |$name| the names of input files can be
c given. All directories default to the current directory,
c whereas all file names are empty, i.e., no input files are
c given.
c
c DOCS	FILENAME	File names
c
c Strings in section |$name| enable the specification of files
c that account for initial conditions or forcing.
c
c |zinit|	Name of file containing initial conditions for water level
c |uvinit|	Name of file containing initial conditions for velocity

	call addfnm('zinit',' ')
        call addfnm('uvinit',' ')

c |wind|	File with wind data. The file may be either
c		formatted or unformatted. For the format of the unformatted
c		file please see the section where the WIN
c		file is discussed.
c		The format of formatted ASCII file
c		is in standard time-series format, with the
c		first column containing the time in seconds and
c		the next two columns containing the wind data.
c		The meaning of the two values depend on the
c		value of the parameter |iwtype| in the |para| section.
c |qflux|	File with heat flux data. This file must be in
c		a special format to account for the various parameters
c		that are needed by the heat flux module to run. Please
c		refer to the information on the file |qflux|.
c |rain|	File with rain data. This file is a standard time series
c		with the time in seconds and the rain values
c		in mm/day. The values may include also evaporation. Therefore,
c		also negative values (for evaporation) are permitted.
c |ice|		File with ice cover. The values range from 0 (no ice cover)
c		to 1 (complete ice cover).

	call addfnm('wind',' ')
        call addfnm('rain',' ')
        call addfnm('qflux',' ')
        call addfnm('ice',' ')

c |surfvel|	File with surface velocities from observation. These
c		data can be used for assimilation into the model.
c |restrt|	Name of the file if a restart is to be performed. The
c		file has to be produced by a previous run
c		with the parameter |idtrst| different
c		from 0. The data record to be used in the file for the
c		restart must be given by time |itrst|.
c |gotmpa|	Name of file containing the parameters for the
c		GOTM turbulence model (iturb = 1).

        call addfnm('surfvel',' ')
	call addfnm('restrt',' ')
	call addfnm('gotmpa',' ')

c |tempin|	Name of file containing initial conditions for temperature
c |saltin|	Name of file containing initial conditions for salinity
c |conzin|	Name of file containing initial conditions for concentration

        call addfnm('tempin',' ')
        call addfnm('saltin',' ')
        call addfnm('conzin',' ')

c |tempobs|	Name of file containing observations for temperature
c |saltobs|	Name of file containing observations for salinity

        call addfnm('tempobs',' ')
        call addfnm('saltobs',' ')

c |temptau|	Name of file containing the time scale for nudging 
c		of temperature
c |salttau|	Name of file containing the time scale for nudging 
c		of salinity

        call addfnm('temptau',' ')
        call addfnm('salttau',' ')

c |bfmini|	Name of file containing initial conditions for bfm
c |offlin|	Name of the file if a offline is to be performed. The
c		file has to be produced by a previous run
c		with the parameter |idtoff| greater than 0.

        call addfnm('bfmini',' ')
	call addfnm('offlin',' ')

c DOCS	END

cc non-documented -> try first	HACK	-> initial conditions

        call addfnm('bioin',' ')
        call addfnm('biosin',' ')
        call addfnm('toxi',' ')
        call addfnm('mercin',' ')	!mercury


c |petsc_zcfg|	Name of file containing the configuration of 
c		PETSc solver for zeta
c |amgx_zcfg|	Name of file containing the configuration of 
c		AmgX solver for zeta

	call addfnm('petsc_zcfg',' ')
	call addfnm('amgx_zcfg',' ')

cc ACQUBC

	call fnm_aquabc_init

cc DOCS	DELWAQ		delwaq model

	call addfnm('voldwq',' ')
	call addfnm('flowdwq',' ')
	call addfnm('areadwq',' ')
	call addfnm('femdwq',' ')
	call addfnm('pntfem',' ')
	call addfnm('ftodwq',' ')

cc DOCS	INTERNAL	internal administration - blank section

	call sctfnm(' ')		!sets default section

	call addfnm('title',' ')

	if( .not. haspar('basnam') ) call addfnm('basnam',' ')
	if( .not. haspar('runnam') ) call addfnm('runnam',' ')
	if( .not. haspar('apnnam') ) call addfnm('apnnam',' ')

	call addfnm('apnfil','apnstd.str')
	call addfnm('memfil','.memory')

	call addfnm('pltfil','plot')

	end

c************************************************************************

        subroutine fnm_aquabc_init

        implicit none

cc for model aquabc (curonian)

        call addfnm('biocon',' ')
        call addfnm('bioscon',' ')
        call addfnm('biolight',' ')
        call addfnm('bioaow',' ')  !ascii output for WC
        call addfnm('bioaos',' ')  !ascii output for BS

        !call addfnm('bioph',' ')
        !call addfnm('biotemp',' ')

        call addfnm('bioload',' ') ! point source loads, not tested yet
        call addfnm('bbs_lev',' ') ! BS depth levels
        call addfnm('settl',' ')   ! settling velocities
        call addfnm('ldbs_m',' ')  ! nutrient load from BS A/D forcing, mud
        call addfnm('ldbs_s',' ')  ! nutrient load from BS A/D forcing, sand
        call addfnm('ldrdx_m',' ') ! nutrient load from WC redox forcing, mud

        end

c************************************************************************

