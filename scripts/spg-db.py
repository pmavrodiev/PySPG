#!/usr/bin/python


import spg.utils as utils
from spg.parameter import EnsembleBuilder, ParameterEnsemble
from spg.master import MasterDB
from spg.cmdline import BaseDBCommandParser
from spg import VAR_PATH, RUN_DIR

import sqlite3 as sql
import sys, optparse
import os, os.path


import fnmatch


class DBCommandParser(BaseDBCommandParser):
    """DB command handler"""

 
    def __init__(self):
        BaseDBCommandParser.__init__(self)
        self.prompt = "| spg-db :::~ "
        self.possible_keys = ['weight', 'repeat',  'status', 'queue']



    def do_init(self, c):
        """build PARAMETERS_NAME|DB_NAME [VAR1=VALUE1[:VAR2=VALUE2]]
        Generates a new database out of a parameters.dat"""
        c = c.split()
        i_arg = c[0]
        try:
            i_arg, db_name = self.translate_name(i_arg)
        except: 
            utils.newline_msg("ERR", "results db '%s' doesn't exist. Can not reinit it" )
            return
#        if self.master_db.result_dbs.has_key( db_name ):
#            utils.newline_msg("WRN", "results db '%s' already registered"%self.shorten_name( db_name ), 2)
#            os.remove(db_name)
#            self.do_remove(i_arg) 

        self.current_param_db = ParameterEnsemble( db_name, init_db = False ) 
        if len(c) >1: self.do_set( ":".join( c[1:] ) )

        parser = EnsembleBuilder( stream = open(i_arg), db_name=db_name  )
        parser.init_db(  )
        parser.fill_status(repeat = self.current_param_db.repeat ) 
        self.master_db.update_result_db( self.current_param_db )
        self.current_param_db.id = self.master_db.cursor.lastrowid
        self.master_db.initialise_result_dbs()


    def do_reinit(self, c):
        """build PARAMETERS_NAME|DB_NAME [VAR1=VALUE1[:VAR2=VALUE2]]
        Generates a new database out of a parameters.dat"""
        c = c.split()
        i_arg = c[0]
        try:
            i_arg, db_name = self.translate_name(i_arg)
        except: 
            utils.newline_msg("ERR", "results db '%s' doesn't exist. Can not reinit it" )
            return
        
        #if self.master_db.result_dbs.has_key( db_name ):
       #     utils.newline_msg("WRN", "results db '%s' already registered"%self.shorten_name( db_name ), 2)
        os.remove(db_name)
        self.do_remove(i_arg) 

        self.current_param_db = ParameterEnsemble( db_name, init_db = False ) 
        if len(c) >1: self.do_set( ":".join( c[1:] ) )

        parser = EnsembleBuilder( stream = open(i_arg), db_name=db_name  )
        parser.init_db(  )
        parser.fill_status(repeat = self.current_param_db.repeat ) 
        self.master_db.update_result_db( self.current_param_db )
        self.current_param_db.id = self.master_db.cursor.lastrowid
        self.master_db.initialise_result_dbs()

    def complete_init(self, text, line, begidx, endidx):    
        completions = fnmatch.filter( os.listdir("."), "results*.sqlite" )
        completions.extend( fnmatch.filter( os.listdir("."), "parameters*.dat" ) )
        if text:
            completions = [ f
                            for f in completions
                            if f.startswith(text)
                            ]
        return completions
        
    def do_register(self,c):
        """registers a given results database into the master database"""
        c = c.split()
        i_arg = c[0]
        i_arg, db_name = self.translate_name(i_arg)
        if self.master_db.result_dbs.has_key( db_name ):
            utils.newline_msg("WRN", "results db '%s' already registered"%self.shorten_name( db_name ), 2)
            return 

        self.current_param_db = ParameterEnsemble( db_name ) 
        if len(c) >1: 
            self.do_set( ":".join( c[1:] ) )

        self.master_db.update_result_db( self.current_param_db )
        self.current_param_db.id = self.master_db.cursor.lastrowid
        self.master_db.initialise_result_dbs()
        print " *--- registered '%s'with id=%d  "%(   self.current_param_db.full_name, self.current_param_db.id )
        

    def do_clean(self, c):
        """cleans the results database by setting the R into N"""
        if self.current_param_db:
            self.current_param_db.execute_query('UPDATE run_status SET status = "N" WHERE status ="R"')


    def do_clean_all(self, c):
        """cleans the results database by setting all into N"""
        if self.current_param_db:
            self.current_param_db.execute_query('UPDATE run_status SET status = "N"  ')

    def do_remove(self, c):
        """removes one results database from the registered ones"""
        if not c:
            ls_res_db = [ self.current_param_db.full_name ]
        else:
            ls_res_db = self.filter_db_list( filter = c )
        if not ls_res_db: return

        self.current_param_db = None
        
        for i in ls_res_db: 
            # print i
            self.master_db.execute_query("DELETE FROM dbs WHERE full_name = ?", i  )
            del self.master_db.result_dbs[i]
        self.master_db.synchronise_master()
 
    def complete_remove(self, text, line, begidx, endidx):
        return self.complete(text)
    
    def do_set(self, c):
        """sets a VAR1=VALUE1[:VAR2=VALUE2]
        sets a value in the currently loaded database """
        if not self.current_param_db: 
            utils.newline_msg("WRN", "current db not set... skipping")
            return 
        
        ret = utils.parse_to_dict(c, allowed_keys=self.possible_keys)
        if not ret: 
            return
        
        for k in ret.keys():
            self.current_param_db.__dict__[k] = ret[k]
            if k == "repeat": continue # repeat is not in the master db (should it be added)
            self.master_db.execute_query( 'UPDATE dbs SET %s= ? WHERE id = ?'%k, ret[k], self.current_param_db.id  )


    def __set_status(self, c, st):
        if not c:
            ls_res_db = [ self.current_param_db.full_name ]
        else:
            ls_res_db = self.filter_db_list( filter = c )
        if not ls_res_db: return
        
        for i in ls_res_db: 
            self.current_param_db.status = st
            self.master_db.execute_query( 'UPDATE dbs SET status= ? WHERE id = ?', st, self.current_param_db.id  )

    def do_stop(self, c):
        """stops the currently loaded registered database"""
        self.__set_status(c, 'S')
                
    def do_start(self, c):
        """starts the currently loaded registered database"""
        self.__set_status(c, 'R')
                
    def do_pause(self, c):
         """pauses the currently loaded registered database"""
         self.__set_status(c, 'P')
                
    

if __name__ == '__main__':
    cmd_line = DBCommandParser()
    if len(sys.argv) == 1:
        cmd_line.cmdloop()
    else:
        cmd_line.onecmd(" ".join(sys.argv[1:]))
        


