"""Test script for module :mod:`sonat.render`"""

import util

from sonat.render import render_html_template

def test_render_html_template_base():
    """Just test that it does not fail"""

    # Embedded mode
    shtmle = render_html_template('base.html')

    # Standalone mode
    shtmls = render_html_template('base.html', standalone=True)

    # Checks
    assert len(shtmls) > len(shtmle)

def test_render_html_template_dict2tree():
    """Just test that it does not fail"""

    # Define content
    content = {
        'Single': 'figure.png',
        'Lists': ['item1',  'item2'],
        'Dict':{
            'key1':'value1',
            'key2':'value2',
            },
        'Sub':{
            'SubDict':{
                'subkey1':'subvalue1'
                }
            }
        }

    # Render
    shtml = render_html_template('dict2tree.html', content=content)

    # Checks
    assert '<img src="figure.png"/>' in shtml
    assert '<span class="d2tItemName">key2</span><br/>' in shtml
    assert """<li>
            <span class="d2tSectionName">Lists</span><br/>""" in shtml
    assert """<li>
                    item1""" in shtml

if __name__=='__main__':
    test_render_html_template_base()
    test_render_html_template_dict2tree()
