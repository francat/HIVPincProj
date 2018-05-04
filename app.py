import os
from flask import Flask, render_template, request

app = Flask(__name__)

APP_ROOT= os.path.dirname(os.path.abspath(__file__))

@app.route("/")
def index():
	return render_template("uploadButton.html")

@app.route("/upload")
def upload():
	target = os.path.join(APP_ROOT, '/files/')
	print(target)

	if not os.path.isdir(target):
		os.mkdir(target)

	for file in request.files.getlist("file"):
		print(file)
		filename = file.filename
		destination = "/".join([target,filename])
		print(destination)
		file.save(destination)
	return render_template(uploadButtonl.html)

if __name__ == '__main__':
	app.run()