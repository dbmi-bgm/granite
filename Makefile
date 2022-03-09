configure:
	pip install poetry

.PHONY: build
build:
	poetry install

update:
	poetry update

test:
	pytest -vv

help:
	@make info

clean:
	rm -rf *.egg-info

info:
	@: $(info Here are some 'make' options:)
	   $(info - Use 'make configure' to configure the repo by installing poetry.)
	   $(info - Use 'make build' to install entry point commands.)
	   $(info - Use 'make update' to update dependencies (and the lock file).)
	   $(info - Use 'make test' to run tests.)
