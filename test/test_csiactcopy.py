#! /usr/bin/env python
# ==========================================================================
# This scripts performs unit tests for the csiactcopy script.
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
import gammalib
import cscripts


# ================================ #
# Test class for csiactcopy script #
# ================================ #
class Test(gammalib.GPythonTestSuite):
    """
    Test class for csiactcopy script.

    This test class makes unit tests for the csiactcopy script by using it
    from the command line and from Python.
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

    # Set test functions
    def set(self):
        """
        Set all test functions.
        """
        # Set test name
        self.name("csiactcopy")

        # Append tests
        self.append(self._test_cmd, "Test csiactcopy on command line")
        self.append(self._test_python, "Test csiactcopy from Python")

        # Return
        return

    # Test csiactcopy on command line
    def _test_cmd(self):
        """
        Test csiactcopy on the command line.
        """
        # Kluge to set the command (installed version has no README file)
        if os.path.isfile("README.md"):
            csiactcopy = "../cscripts/csiactcopy.py"
        else:
            csiactcopy = "csiactcopy"

        # Setup csiactcopy command
        cmd = csiactcopy+' remote_master="iactdata/master.json"'+ \
                         ' prodname="unit-test"'+ \
                         ' outpath="iactdata_cmd1"'+ \
                         ' logfile="csiactcopy_cmd1.log" chatter=1'

        # Execute csiactcopy, make sure we catch any exception
        try:
            rc = os.system(cmd+" >/dev/null 2>&1")
        except:
            pass

        # Check if execution was successful
        self.test_assert(rc == 0,
                         "Successful csiactcopy execution on command line")

        # Check copy
        self._check_copy("iactdata_cmd1")

        # Setup csiactcopy command
        cmd = csiactcopy+' remote_master="master_that_does_not_exist"'+ \
                         ' prodname="unit-test"'+ \
                         ' outpath="iactdata_test"'+ \
                         ' logfile="csiactcopy_cmd2.log"'

        # Execute csiactcopy, make sure we catch any exception
        try:
            rc = os.system(cmd+" >/dev/null 2>&1")
        except:
            pass

        # Check if execution failed
        self.test_assert(rc != 0,
                         "Failure of csiactcopy execution on command line")

        # Return
        return

    # Test csiactcopy from Python
    def _test_python(self):
        """
        Test csiactcopy from Python.
        """
        # Set-up csiactcopy
        findobs = cscripts.csiactcopy()
        findobs["remote_master"] = "iactdata/master.json"
        findobs["prodname"]      = "unit-test"
        findobs["outpath"]       = "iactdata_py1"
        findobs["logfile"]       = "csiactcopy_py1.log"
        findobs["chatter"]       = 2

        # Run csiactcopy script and save run list
        findobs.logFileOpen()   # Make sure we get a log file
        findobs.run()
        findobs.save()

        # Check copy
        self._check_copy("iactdata_py1")

        # Return
        return

    # Check copy
    def _check_copy(self, pathname):
        """
        Check copy.
        """
        # Set file names
        hdu_index_name = gammalib.GFilename(pathname+"/hdu-index.fits")
        obs_index_name = gammalib.GFilename(pathname+"/obs-index.fits")
        master_name    = gammalib.GFilename(pathname+"/master.json")

        # Check for existence of files
        self.test_assert(hdu_index_name.exists(),
                         'Check is file "hdu-index.fits" exists.')
        self.test_assert(obs_index_name.exists(),
                         'Check is file "obs-index.fits" exists.')
        self.test_assert(master_name.exists(),
                         'Check is file "master.json" exists.')

        # Check if index files are FITS file
        self.test_assert(hdu_index_name.is_fits(),
                         'Check is file "hdu-index.fits" is FITS file.')
        self.test_assert(obs_index_name.is_fits(),
                         'Check is file "obs-index.fits" is FITS file.')
        
        # Return
        return
