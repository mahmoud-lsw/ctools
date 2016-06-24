#! /usr/bin/env python
# ==========================================================================
# Performs unit tests for the ctools example scripts
#
# Copyright (C) 2016 Juergen Knoedlseder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==========================================================================
import os
import sys
import gammalib
import ctools
import cscripts


# =============================================== #
# Test class for example executables unit testing #
# =============================================== #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for GammaLib example executables
    """
    # Constructor
    def __init__(self):
        """
        Constructor.
        """
        # Call base class constructor
        gammalib.GPythonTestSuite.__init__(self)

        # Return
        return

    # Execute Python script
    def _execute_python(self, name, args='', log=''):
        """
        Execute Python script
        
        Parameters
        ----------
        name : str
            Python script name with .py extension
        args : str, optional
            String with arguments to be passed to script
        log : str, optional
            String with name to be used for logging
        """
        # Setup command
        cmd = '../examples/' + name + '.py'
        if len(args) > 0:
            cmd += ' ' + args

        # Set logname
        if len(log) > 0:
            logname = log
        else:
            logname = name

        # Execute Python script, make sure we catch any exception
        try:
            rc = os.system(cmd+' > example_'+logname+'.log 2>&1')
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0, 'Check "'+name+'" script from command line')

        # Return
        return

    # Set test functions
    def set(self):
        """
        Set all test functions
        """
        # Set test name
        self.name('examples')

        # Append tests
        self.append(self.test_generate_prod3_irfs, 'Test generate_prod3_irfs.py')
        self.append(self.test_make_pointings, 'Test make_pointings.py')
        self.append(self.test_make_spectrum, 'Test make_spectrum.py')
        self.append(self.test_make_ts_distributions, 'Test make_ts_distributions.py')
        self.append(self.test_pipeline_binned_disk, 'Test pipeline_binned_disk.py')
        self.append(self.test_pipeline_binned_mem, 'Test pipeline_binned_mem.py')
        self.append(self.test_pipeline_stacked_disk, 'Test pipeline_stacked_disk.py')
        self.append(self.test_pipeline_stacked_mem, 'Test pipeline_stacked_mem.py')
        self.append(self.test_pipeline_unbinned_disk, 'Test pipeline_unbinned_disk.py')
        self.append(self.test_pipeline_unbinned_mem, 'Test pipeline_unbinned_mem.py')
        self.append(self.test_show_butterfly, 'Test show_butterfly.py')
        self.append(self.test_show_model, 'Test show_model.py')
        self.append(self.test_show_pha, 'Test show_pha.py')
        self.append(self.test_show_pull_evolution, 'Test show_pull_evolution.py')
        self.append(self.test_show_pull_histogram, 'Test show_pull_histogram.py')

        # Return
        return

    # Test make_pointings
    def test_make_pointings(self):
        """
        Test make_pointings
        """
        # Set script name
        script = 'make_pointings'

        # Execute script
        self._execute_python(script, args='-h', log='make_pointings_h')
        self._execute_python(script, args='gps', log='make_pointings_gps')
        self._execute_python(script, args='gps3', log='make_pointings_gps3')
        self._execute_python(script, args='extgal', log='make_pointings_extgal')
        self._execute_python(script, args='gc', log='make_pointings_gc')
        self._execute_python(script, args='lmc', log='make_pointings_lmc')

        #TODO: Do any testing

        # Return
        return

    # Test generate_prod3_irfs
    def test_generate_prod3_irfs(self):
        """
        Test generate_prod3_irfs
        """
        # Execute script
        self._execute_python('generate_prod3_irfs')

        #TODO: Do any testing

        # Return
        return

    # Test make_spectrum
    def test_make_spectrum(self):
        """
        Test make_spectrum
        """
        # Execute script
        self._execute_python('make_spectrum', args='eps')

        #TODO: Do any testing

        # Return
        return

    # Test make_ts_distributions
    def test_make_ts_distributions(self):
        """
        Test make_ts_distributions
        """
        # Execute script
        self._execute_python('make_ts_distributions', args='-n 2 -e 0 -d 1800')

        #TODO: Do any testing

        # Return
        return

    # Test pipeline_binned_disk
    def test_pipeline_binned_disk(self):
        """
        Test pipeline_binned_disk
        """
        # Execute script
        self._execute_python('pipeline_binned_disk')

        #TODO: Do any testing

        # Return
        return

    # Test pipeline_binned_mem
    def test_pipeline_binned_mem(self):
        """
        Test pipeline_binned_mem
        """
        # Execute script
        self._execute_python('pipeline_binned_mem')

        #TODO: Do any testing

        # Return
        return

    # Test pipeline_stacked_disk
    def test_pipeline_stacked_disk(self):
        """
        Test pipeline_stacked_disk
        """
        # Execute script
        self._execute_python('pipeline_stacked_disk')

        #TODO: Do any testing

        # Return
        return

    # Test pipeline_stacked_mem
    def test_pipeline_stacked_mem(self):
        """
        Test pipeline_stacked_mem
        """
        # Execute script
        self._execute_python('pipeline_stacked_mem')

        #TODO: Do any testing

        # Return
        return

    # Test pipeline_unbinned_disk
    def test_pipeline_unbinned_disk(self):
        """
        Test pipeline_unbinned_disk
        """
        # Execute script
        self._execute_python('pipeline_unbinned_disk')

        #TODO: Do any testing

        # Return
        return

    # Test pipeline_unbinned_mem
    def test_pipeline_unbinned_mem(self):
        """
        Test pipeline_unbinned_mem
        """
        # Execute script
        self._execute_python('pipeline_unbinned_mem')

        #TODO: Do any testing

        # Return
        return

    # Test show_butterfly
    def test_show_butterfly(self):
        """
        Test show_butterfly
        """
        # Execute script
        self._execute_python('show_butterfly',
                             args='data/butterfly.txt example_butterfly.eps')

        #TODO: Do any testing

        # Return
        return

    # Test show_model
    def test_show_model(self):
        """
        Test show_model
        """
        # Execute script
        self._execute_python('show_model',
                             args='data/crab.xml Crab example_model.eps')

        #TODO: Do any testing

        # Return
        return

    # Test show_pha
    def test_show_pha(self):
        """
        Test show_pha
        """
        # Execute script
        self._execute_python('show_pha', args='data/pha.fits example_pha.eps')

        #TODO: Do any testing

        # Return
        return

    # Test show_pull_evolution
    def test_show_pull_evolution(self):
        """
        Test show_pull_evolution
        """
        # Execute script
        self._execute_python('show_pull_evolution',
                             args='data/pull.dat Pull_Crab_Prefactor '
                                  'example_pull_evolution.eps')

        #TODO: Do any testing

        # Return
        return

    # Test show_pull_histogram
    def test_show_pull_histogram(self):
        """
        Test show_pull_histogram
        """
        # Execute script
        self._execute_python('show_pull_histogram',
                             args='data/pull.dat Pull_Crab_Prefactor 50 '
                                  'example_pull_histogram.eps')

        #TODO: Do any testing

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Set calibration database
    os.environ['CALDB'] = '../caldb'

    # Set PFILES environment variable
    try:
        os.mkdir('pfiles')
    except:
        pass
    os.environ['PFILES'] = 'pfiles'

    # Copy over pfiles
    os.system('cp -r ../src/*/*.par pfiles/')
    os.system('cp -r ../cscripts/*.par pfiles/')

    # Allocate test suites
    suites = gammalib.GTestSuites('Examples testing')

    # Allocate test suite, setup tests and append them to the container
    suite = Test()
    suite.set()
    suites.append(suite)

    # Run test suite
    success = suites.run()

    # Save test results
    suites.save('reports/examples.xml')

    # Set return code
    if success:
        rc = 0
    else:
        rc = 1

    # Exit with return code
    sys.exit(rc)
