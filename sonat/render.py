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
        <script type="text/javascript" src="CollapsibleLists.compressed.js"></script>
        <script type="text/javascript">

          runOnLoad(function(){ CollapsibleLists.apply(); });

        </script>
        <link rel="stylesheet" href="sonat.css" type="text/css">
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
    <ul class="treeView collapsibleList">
    {% for name0, content0 in content.iteritems() %}
        {% if content0  %}
        <li>
            <span class="d2tSectionName">{{ name0 }}</span><br/>

            {% if content0 is string %}
            {{ content0|checkimg }}
            {% elif content0 is mapping %}
            <ul class="xxxcollapsibleList">
            {% for name1, content1 in content0.iteritems() %}
                {% if content1 %}
                <li>
                    <span class="d2tItemName">{{ name1 }}</span><br/>

                    {% if content1 is string %}
                    {{ content1|checkimg }}
                    {% elif content1 is mapping %}
                    <ul>
                    {% for name2, content2 in content1.iteritems() %}
                        {% if content2 %}
                        <li>
                            <span class="d2tItemName">{{ name2 }}</span><br/>

                            {% if content2 is string %}
                            {{ content2|checkimg }}
                            {% elif content2 is mapping %}
                            <ul>
                            {% for name3, content3 in content2.iteritems() %}
                                {% if content3 %}
                                <li>
                                    <span class="d2tItemName">{{ name3 }}</span><br/>

                                    {{ content3|checkimg }}

                                </li>
                                {% endif %}
                            {% endfor %}
                            </ul>
                            {% elif content2 is sequence %}
                            <ul>
                            {% for item in content2 %}
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
                    {% elif content1 is sequence %}
                    <ul>
                    {% for item in content1 %}
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
            {% elif content0 is sequence %}
            <ul class="collapsibleList">
            {% for item in content0 %}
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

#: Html dependecies
HTML_DEPS = ["runOnLoad.js", "CollapsibleLists.compressed.js", "sonat.css",
    "button-closed.png",
    "button-open.png",
    "button.png",
    "list-item-contents.png",
    "list-item-last-open.png",
    "list-item-last.png",
    "list-item-open.png",
    "list-item-root.png",
    "list-item.png",
]

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
    for bfile in HTML_DEPS:
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

    # Check material
    copy_html_material(os.path.dirname(htmlfile))



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

