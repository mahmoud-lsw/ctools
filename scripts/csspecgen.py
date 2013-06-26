#! /usr/bin/env python
import gammalib

class csspecgen(gammalib.GApplication):
    """
    TODO: document csspecgen tool.
    """
    def __init__(self, *argv):
        self.name = "csspecgen"
        self.version = "0.1.0"
        
        # Initialise some members
        self.m_obs = None
        
        # Make sure that parfile exists
        file = self.parfile()

        # Initialise application
        if len(argv) == 0:
            gammalib.GApplication.__init__(self, self.name, self.version)
        elif len(argv) == 1:
            gammalib.GApplication.__init__(self, self.name, self.version, *argv)
        else:
            raise TypeError("Invalid number of arguments given.")

        # Set logger properties
        self.log_header()
        self.log.date(True)
    
    def __del__(self):
        #  Write separator into logger
        if self.logTerse():
            self.log("\n")
        
    def get_parameters(self):
        """
        Get parameters from parfile and setup the observation.
        """
        # Get parameters
        self.m_evfile = self["evfile"].filename()
        self.m_onregfile = self["onregfile"].filename()
        self.m_exregfile = self["exregfile"].filename()
        self.m_outphaprefix = self["outphaprefix"].filename()
        self.m_outregfile = self["outregfile"].filename()

        self.m_obs = gammalib.GCTAObservations()
        # Try first to open as FITS file
        try:
            # Load event list in CTA observation
            obs = gammalib.GCTAObservation()
            obs.load_unbinned(m_evfile);
            m_obs.append(obs);
        except gammalib.GException.fits_open_error:
            # Load observations from XML file
            m_obs.load(m_evfile);

        self.m_emin = self["emin"].real()
        self.m_emax = self["emax"].real()
        self.m_enumbins = self["enumbins"].integer()

        # Set some fixed parameters
        self.m_log = False  # Logging in client tools
        self.m_debug = False  # Debugging in client tools
        
    def execute(self):
        self.run()        
        self.m_onoff_obs.save_regions(self.m_outregfile)
        self.m_onoff_obs.save_pha(self.m_outphaprefix)

    def run(self):
        # Switch screen logging on in debug mode
        if self.logDebug():
            self.log.cout(True)

        # Get parameters
        self.get_parameters()
        
        #  Write input parameters into logger
        if self.logTerse():
            self.log_parameters()
            self.log("\n")
        
        # Write observation into logger
        if self.logTerse():
            self.log("\n")
            self.log.header1("Observation")
            self.log(str(self.m_obs))
            self.log("\n")

        # Write header
        if self.logTerse():
            self.log("\n")
            self.log.header1("Input")


        self.m_on_region = gammalib.GRegions.load(self.m_onregfile)
        self.log(str(self.m_on_region))

        self.m_exclusion_regions = gammalib.GRegions.load(self.m_exregfile)

        # m_obs has been set in get_parameters already

        # Setup energy range covered by data
        emin = gammlib.GEnergy(m_emin, "TeV");
        emax = gammalib.GEnergy(m_emax, "TeV");
        self.m_ebds = gammalib.GEbounds(m_enumbins, emin, emax);
        
        self.m_onoff_obs = gammalib.GCTAOnOffObservations(self.m_obs, self.m_ebds)
        
        for spec in self.m_onoff_obs.spectra:
            # spec is a GCTAOnOffSpectrum
            spec.set_on_region(self.m_on_region)
            spec.set_exclusion_regions(self.m_exclusion_regions)
            spec.compute_off_regions()
            spec.fill_spectra() # fills n_on, n_off

if __name__ == '__main__':
    # Create instance of application
    app = csspecgen(sys.argv)
    
    # Open logfile
    app.logFileOpen()
    
    # Execute application
    app.execute()
