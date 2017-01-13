"""Rendering of results"""
# -*- conding: utf8 -*-

import os
import shutil
import re
from jinja2 import Template, Environment, DictLoader
from jinja2.filters import FILTERS as JINJA_FILTERS
from jinja2.utils import soft_unicode


TEMPLATE_HTML_BASE = """
<!DOCTYPE html>
<html>

    <head>
        <script type="text/javascript" src="runOnLoad.js"></script>
        <script type="text/javascript" src="CollapsibleLists.js"></script>
        <script type="text/javascript">

          runOnLoad(function(){ CollapsibleLists.apply(); });

        </script>
        <link rel="stylesheet" href="render.css" type="text/css">
        <title>SONAT - {{ title }}</title>
    </head>

    <body>
    {% if standalone %} <h1>{{ title }}</h1> {% endif %}
    {% block content %}
    {% endblock %}
    {% if standalone %} <div id="footer">Created with SONAT</div> {% endif %}
    </body>

</html>
"""

TEMPLATE_HTML_DICT2TREE = """
{% extends "base.html" %}

{% block content %}
    <ul>
    {% for section_name, subcontent in content.iteritems() %}
        {% if subcontent  %}
        <li>
            <span class="d2tSectionName">{{ section_name }}</span><br/>
            {% if subcontent is string %}
            {{ subcontent|checkimg }}
            {% elif subcontent is mapping %}
            <ul>
            {% for item_name, item_content in subcontent.iteritems() %}
                {% if item_content %}
                <li>
                    <span class="d2tItemName">{{ item_name }}</span><br/>
                    {{ item_content|checkimg }}
                </li>
                {% endif %}
            {% endfor %}
            </ul>
            {% elif subcontent is sequence %}
            <ul>
            {% for item in subcontent %}
                {% if item %}
                <li>
                    {{ item|checkimg }}
                </li>
                {% endif %}
            {% endfor %}
            </ul>

            {% endif %}
        </li>
        {% endif %}
    {% endfor %}
  </ul>

{% endblock %}
"""


#: Mapping for the jinja loader
JINJA_LOADER_MAPPING = {}

#: Jinja loader
JINJA_LOADER = DictLoader(JINJA_LOADER_MAPPING)

#: Jinja environment
JINJA_ENV = Environment(
    trim_blocks=True,
    lstrip_blocks=True,
    loader=JINJA_LOADER)

def copy_html_material(workdir=None):
    """Copy all the needed material for html rendering to the working dir"""

    # Reference directories
    thisdir = os.path.dirname(__file__)
    if workdir is None:
        workdir = os.getcwd()

    # Copies
    for bfile in ["runOnLoad.js", "CollapsibleLists.js", "sonat.css"]:
        sfile = os.path.join(thisdir, bfile)
        dfile = os.path.join(workdir, bfile)
        if not os.path.exists(dfile):
            shutil.copy(sfile, dfile)


def register_html_template(name, content):
    """Register a jinja template into the :attr:`TEMPLATES` variable"""
    JINJA_LOADER_MAPPING[name] = content
#    HTML_TEMPLATES[name] = Template(content)


def render_html_template(name, **kwargs):
    """Render a jinja template"""
#    kwargs.update(HTML_TEMPLATES)
    kwargs.setdefault('standalone', False)
#    return HTML_TEMPLATES[name].render(**kwargs)
    return JINJA_ENV.get_template(name).render(**kwargs)

def render_and_export_html_template(name, htmlfile, **kwargs):
    """Render a jinja template and export it to an html file"""
    # Export to an html string
    shtml = render_html_template(name, **kwargs)

    # Write to file
    f = open(htmlfile, 'w')
    f.write(shtml.encode('utf8'))
    f.close


RE_IMG_MATCH = re.compile(r'.*(png|jpg|jpeg)$', re.I).match
def do_checkimg(obj):
    sobj = soft_unicode(obj)
    if RE_IMG_MATCH(sobj):
        return '<img src="{}"/>'.format(sobj)
    return obj
JINJA_ENV.filters['checkimg'] = do_checkimg


# 0HTML_TEMPLATES = {}
register_html_template('base.html', TEMPLATE_HTML_BASE)
register_html_template('dict2tree.html', TEMPLATE_HTML_DICT2TREE)

