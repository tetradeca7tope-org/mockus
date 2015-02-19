% unit test for project to Rectangle.

pt = rand(5, 1), projectToRectangle(pt),

bounds = [1 5; 3 10; -9 -3];
pt = [2; 11; -10]; [bounds pt projectToRectangle(pt, bounds)],
pt = [2; 6; -10]; [bounds pt projectToRectangle(pt, bounds)],

