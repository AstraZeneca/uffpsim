build-dev:
	docker build -t uffpsim .

run-dev:
	docker run -v `pwd`:/workspace --rm -it uffpsim bash
