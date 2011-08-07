#!/usr/bin/python


from spg.parameter import ParameterExecutor
from spg import BINARY_PATH


import os, os.path
import optparse


if __name__ == "__main__":

    parser = optparse.OptionParser(usage = "usage: %prog [options] project_id1 project_id2 project_id3... ")
    parser.add_option("--timeout", type="int", action='store', dest="timeout",
                            default = 60 , help = "timeout for database connection" )

    parser.add_option("--tree", action='store_true', dest="tree",
                       help = "whether to create a directory tree with the key-value pairs" )

    parser.add_option("-d","--directory-var", action='store', type = "string", dest="directory_vars",
                       default = False, help = "which variables to store as directories, only if tree" )

    options, args = parser.parse_args()
    
    if len(args) == 0:
        args = ["results.sqlite"]
    
    for i_arg in args:
      if ".sqlite" not in i_arg:
          db_name = i_arg.replace("parameters","").replace(".dat","")
          db_name = "results%s.sqlite"%db_name
      else:
          db_name = i_arg
      full_name = os.path.realpath(db_name)
#      path, out = os.path.split(full_name)
      executor = ParameterExecutor( full_name )
      
      executor.generate_tree( options.directory_vars )

      for values in executor:
          executor.launch_process()
          
          #      parser.init_db()
#          parser.fill_status(repeat = options.repeat )

