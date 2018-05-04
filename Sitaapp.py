import os
from flask import Flask, render_template, request


app = Flask(__name__)

APP_ROOT = os.path.dirname(os.path.abspath(__file__))


@app.route("/upload", methods = ['POST'])
def upload():
	return render_template('uploadButton.html')
	target = os.path.join(APP_ROOT, "files/")
	print(target)


	if not os.path.isdir(target):
		os.mkdir(target)

	for file in request.files.getlist("file"):
		print(file)
		filename = file.filename
		destination = "/".join([target, filename])
		print(destination)
		file.save(destination)
	return render_template("complete.html", filename=filename, destination=destination)
			#                   ^^^ should render algorithmSelectionPage.html

#below: onClick with "run" button on algorithmSelectionPage.html
@app.route("/run, <filename>, <destination>", methods = ['POST'])
def run(filename, destination):
	with open(os.path.join(destination), 'w+') as txt:
		content = txt.read()
		with open(os.path.join(destination), 'w+') as outputFile:
			#line 39: will be replaced with feature making and ML algorithm
			outputFile.write(content.lower())
			outputFile.save(destination)
	target = os.path.join(APP_ROOT, "files_ML/")
	return render_template("complete.html")

if __name__ == "__main__":
	app.run(debug= True)
