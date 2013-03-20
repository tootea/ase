.. _devel:

===============
ASE development
===============

As a developer, you should subscribe to all ASE related
:ref:`mailing_lists`.


Development topics
==================

.. toctree::

   releasenotes
   contribute
   versioncontrol
   python_codingstandard
   writing_documentation_ase
   calculators
   making_movies
   newrelease
   tests
   buildbot
   translate
   todo
   py3k
   proposals/proposals


.. _devel_passwd_create:

Creating an encrypted password for SVN access
=============================================

Use this command::

  htpasswd -nm <your-desired-user-name>

and type a good password twice.  The encrypted password will be
printed on the screen.

If you don't have the ``htpasswd`` command, then use Python:

>>> import crypt
>>> passwd = '<your-password>'
>>> print crypt.crypt(passwd, passwd)
