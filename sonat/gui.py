#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright or © or Copr. Actimar/IFREMER (2010-2017)
#
# This software is a computer program whose purpose is to provide
# utilities for handling oceanographic and atmospheric data,
# with the ultimate goal of validating the MARS model from IFREMER.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#

import collections, pprint, os, sys, traceback

from PyQt4 import QtCore, QtGui, QtWebKit

from vacumm.misc.log import Logger
from vacumm.misc.cfgui.application import Application, Plugin
from vacumm.misc.cfgui import QtObject
from vacumm.misc.cfgui.utils.ui import info_dialog, warning_dialog, error_dialog, confirm_dialog

from sonat.cui import (
    ens_gen_pseudo_from_cfg,
    ens_plot_diags_from_cfg,
    obs_plot_from_cfg,
    arm_analysis_from_cfg,
    arm_sa_from_cfg,
)

# This is at least required to setup sonat config validators
import sonat.config


def find_config(cfg, *path):
    return find_config(cfg[path[0]], *path[1:]) if len(path)>1 else cfg[path[0]]


class SonatPlugin(Plugin):
    
    def __init__(self, *args, **kwargs):
        Plugin.__init__(self, *args, **kwargs)
        
        self.tools = collections.OrderedDict((
            ('ens_gen_pseudo', dict(
                title='Generate a pseudo ensemble',
                func=ens_gen_pseudo_from_cfg,
                type='generate_and_check_file',
                check_file_cfgpath=('ens', 'ncensfile')
            )),
            ('ens_plot_diags', dict(
                title='Plot ensemble diagnostics',
                func=ens_plot_diags_from_cfg,
                type='generate_and_open_html',
                open_html_cfgpath=('ens', 'htmlfile')
            )),
            ('obs_plot', dict(
                title='Plot observation platform locations',
                func=obs_plot_from_cfg,
                type='generate_and_open_html',
                open_html_cfgpath=('obs', 'htmlfile')
            )),
            ('arm_analysis', dict(
                title='Run an ARM analysis',
                func=arm_analysis_from_cfg,
                type='generate_and_open_html',
                open_html_cfgpath=('arm', 'analysis', 'htmlfile')
            )),
            ('arm_sa', dict(
                title='Run ARM sensitivity analyses',
                func=arm_sa_from_cfg,
                type='generate_and_open_html',
                open_html_cfgpath=('arm', 'sa', 'htmlfile')
            )),
            ('sonat_logs', dict(
                title='Show log file',
                func=self.show_log
            )),
        ))
        
        self.tab_widgets = set()
    
    def enable(self):
        Plugin.enable(self)
        
        self.main_window = self.application.main_controller.main_window
        
        self.application.sessions_controller.sessions_dialog.line_specification.setEnabled(False)
        self.application.sessions_controller.sessions_dialog.button_specification.setEnabled(False)
        
        self.menu = QtGui.QMenu(self.main_window.menubar)
        self.menu.setObjectName('menu_sonat')
        self.menu.setTitle('Sonat')
        self.main_window.menubar.insertAction(self.main_window.menu_help.menuAction(), self.menu.menuAction())
        
        for name,tool in self.tools.items():
        
            action = QtGui.QAction(self.main_window)
            action.setIcon(QtGui.QIcon.fromTheme(tool.get('icon', 'application-x-executable')))
            action.setText(tool['title'])
            
            def connect_start_tool(tool):
                def start_tool():
                    self.start_tool(tool)
                action.triggered.connect(start_tool)
            connect_start_tool(tool)
            
            self.menu.addAction(action)
        
        action = QtGui.QAction(self.main_window)
        action.setIcon(QtGui.QIcon.fromTheme(tool.get('icon', 'help-contents')))
        action.setText('User\'s guide')
        action.triggered.connect(self.on_menu_help)
        self.main_window.menu_help.addAction(action)
        
        self.main_window.tabs.tabCloseRequested.connect(self.on_tab_close_requested)
    
    def disable(self):
        Plugin.disable(self)
        # TODO
    
    def session_created(self, session):
        session.specification_file = sonat.config.SONAT_INIFILE
        self.info('Session created:\n%s', session.to_xml_str())
    
    def on_menu_help(self):
        self.add_html_tab('Help', 'http://relay.actimar.fr/~raynaud/sonat')
    
    def on_tab_close_requested(self, index):
        self.info('on_tab_close_requested')
        widget = self.main_window.tabs.widget(index)
        if widget in self.tab_widgets:
            self.main_window.tabs.removeTab(index)
            self.tab_widgets.remove(widget)
            widget.deleteLater()
    
    def add_html_tab(self, title, url):
        self.info('add_html_tab url: %r', url)
        tab  = QtGui.QWidget()
        layout = QtGui.QVBoxLayout(tab)
        tab_index = self.main_window.tabs.addTab(tab, title + ' - ' + os.path.basename(url))
        self.tab_widgets.add(tab)
        self.main_window.tabs.setCurrentIndex(tab_index)
        webview = QtWebKit.QWebView(tab)
        layout.addWidget(webview)
        if not url.startswith('http://') and not url.startswith('https://'):
            url = 'file://'+url
        webview.load(QtCore.QUrl(url))
    
    def get_configuration(self):
        
        #return self.application.main_controller.get_view_configuration()
        
        from sonat.config import load_cfg
        
        #return load_cfg(self.application.main_controller.session.configuration_file)
        
        return load_cfg(self.application.main_controller.get_view_configuration())
    
    def start_tool(self, tool):
        
        name = self.tools.keys()[self.tools.values().index(tool)]
        
        cfg = self.get_configuration()
        
        self.info('start tool %r %r', name, tool)
        self.debug('cfg: %s', pprint.pformat(cfg.dict()))
        
        if 'type' not in tool:
            tool['func']()
            return
        
        if 'thread' in tool:
            warning_dialog('This tool is already running')
            return
        tool['thread'] = RunThread(name, tool, cfg)
        
        progress = QtGui.QProgressDialog()
        progress.setWindowTitle(tool['title'])
        progress.setLabelText(tool['title'])
        progress.setMinimum(0)
        progress.setMaximum(0)
        progress.setCancelButton(None)
        
        if tool['type'] == 'generate_and_check_file':
            
            check_file = find_config(cfg.dict(), *tool['check_file_cfgpath'])
            check_file = os.path.join(cfg['session']['workdir'], check_file)
            
            if os.path.isfile(check_file):
                if not confirm_dialog('File %r already exists, regenerate file ?'%(check_file,), title=tool['title']):
                    return
                os.remove(check_file)
            
            progress.setLabelText('Generating %s'%(check_file,))
            
            def finished():
                if not os.path.isfile(check_file):
                    msg = 'Error while generating file %r'%check_file
                    self.error(msg)
                    error_dialog(msg, title=tool['title'])
                else:
                    msg = 'File %r generated'%check_file
                    self.info(msg)
                    info_dialog(msg, title=tool['title'])
        
        elif tool['type'] == 'generate_and_open_html':
            
            html_file = find_config(cfg.dict(), *tool['open_html_cfgpath'])
            html_file = os.path.join(cfg['session']['workdir'], html_file)
            
            if os.path.isfile(html_file):
                if not confirm_dialog('File %r already exists, regenerate file ?'%(html_file,), title=tool['title']):
                    self.add_html_tab(tool['title'], html_file)
                    return
                os.remove(html_file)
            
            progress.setLabelText('Generating %s'%(html_file,))
            
            def finished():
                if not os.path.isfile(html_file):
                    msg = 'Error while generating file %r'%html_file
                    self.error(msg)
                    error_dialog(msg, title=tool['title'])
                else:
                    msg = 'File %r generated'%html_file
                    self.info(msg)
                    self.add_html_tab(tool['title'], html_file)
                    info_dialog(msg, title=tool['title'])
        
        def wrap_finished(self, name, tool, finished):
            def wrapped_finished(res, exc, tb):
                self.verbose('tool thread finised %r', name)
                progress.hide()
                thread = tool['thread']
                if exc:
                    detail = tb
                    error_dialog('Error running tool %r: %s'%(tool['title'], exc), title=tool['title'], detail=detail)
                else:
                    finished()
                thread.deleteLater()
                tool.pop('thread')
            return wrapped_finished
        tool['thread'].finished.connect(wrap_finished(self, name, tool, finished))
        
        progress.open()
        self.verbose('start tool thread %r', name)
        tool['thread'].start()
    
    def show_log(self):
        title = self.tools['sonat_logs']['title']
        cfg = self.get_configuration()
        log_file = cfg['logger']['file']
        #log_file = os.path.join(cfg['session']['workdir'], log_file)
        if not os.path.exists(log_file):
            error_dialog('Log file %r does not exists'%(log_file,), title=title)
            return
        
        tab  = QtGui.QWidget()
        layout = QtGui.QVBoxLayout(tab)
        tab_index = self.main_window.tabs.addTab(tab, title + ' - ' + os.path.basename(log_file))
        self.tab_widgets.add(tab)
        self.main_window.tabs.setCurrentIndex(tab_index)
        
        sublayout = QtGui.QVBoxLayout(tab)
        layout.addLayout(sublayout)
        
        refresh_button = QtGui.QPushButton(tab)
        refresh_button.setText('Refresh')
        refresh_button.setIcon(QtGui.QIcon.fromTheme('view-refresh'))
        refresh_button.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum))
        sublayout.addWidget(refresh_button)
        
        auto_refresh_checkbox = QtGui.QCheckBox(tab)
        auto_refresh_checkbox.setText('Auto refresh')
        auto_refresh_checkbox.setCheckState(True)
        sublayout.addWidget(auto_refresh_checkbox)
        auto_refresh_checkbox.setTristate(False)
        
        textedit = QtGui.QTextEdit(tab)
        textedit.setReadOnly(True)
        layout.addWidget(textedit)
        
        def refresh():
            with file(log_file, 'rU') as f:
                textedit.setText(f.read())
                textedit.moveCursor(QtGui.QTextCursor.End)
        refresh()
        refresh_button.clicked.connect(refresh)
        
        def setup_auto_refresh():
            data = dict(last_stat= os.stat(log_file))
            def auto_refresh():
                if not auto_refresh_checkbox.isChecked():
                    return
                stat = os.stat(log_file)
                if stat.st_mtime != data['last_stat'].st_mtime:
                    refresh()
                data['last_stat'] = stat
            
            return auto_refresh
        timer = QtCore.QTimer(tab)
        timer.timeout.connect(setup_auto_refresh())
        timer.start(1500)
        


class RunResult(QtCore.QObject):
    def __init__(self, exc=None):
        self.exc = exc

class RunThread(QtObject, QtCore.QThread):
    
    finished = QtCore.pyqtSignal(object, object, object)
    
    def __init__(self, name, tool, cfg):
        QtObject.__init__(self)
        QtCore.QThread.__init__(self)
        self.name = name
        self.tool = tool
        self.cfg = cfg
    
    def run(self):
        try:
            self.info('run tool %r %r', self.name, self.tool)
            res = self.tool['func'](self.cfg)
            self.info('tool %r returned %r', self.name, res)
        except Exception, exc:
            self.exception('tool %r execution failed: %s', self.name, exc)
            # Traceback must be generated from this thread
            self.finished.emit(None, '%s'%exc, traceback.format_exc(exc))
        else:
            self.finished.emit(res, None, None)


def populate_argparser(parser):
    Application.populate_argparser(parser)
    parser.set_defaults(configuration_directory=os.path.join(os.path.expanduser('~'), '.sonat', 'gui'))


def run_from_args(parser, args, cfg):
    Logger.apply_class_argparser_options(args)
    app = Application()
    app.register_plugin(SonatPlugin)
    return app.run(args)
