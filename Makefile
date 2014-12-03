.PHONY: all clean hooks

all:
	cd g2g; make;
	cd lioamber; make;
	cd liosolo; make;

clean:
	cd g2g; make clean;
	cd lioamber; make clean;
	cd liosolo; make clean;

hooks:
	cp hooks/pre-push.py.hook .git/hooks/pre-push
