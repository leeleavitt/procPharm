# install arangoDB community from website. Asks for a super password.
# start arangoShell and enter the super password.
#pip install pyarango

from pyArango.connection import *
conn = Connection(username="root", password="roselab2")

db = conn["test"]