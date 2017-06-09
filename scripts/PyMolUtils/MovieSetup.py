
from pymol import cmd

cmd.reinitialize();
cmd.run('/Users/joy/src/elfin/src/PyMolUtils/LineUtils.py');
cmd.load('/Users/joy/src/elfin/src/Python/Greedy.py');

print '\n\n\n\n\n\n\n'
print 'Loaded'

# For quick testing
# main(['/Users/joy/src/elfin/bm/l10/6vjex8d.json'])
# main(['/Users/joy/src/elfin/bm/l10/kmb0yfh.json'])
# main(['/Users/joy/src/elfin/bm/fun/M.csv'])
main(['/Users/joy/src/elfin/bm/fun/S.csv', 's=2.50'])