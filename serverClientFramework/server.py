from flask import Flask, request, jsonify, render_template, session, redirect
from flask_wtf import CSRFProtect

csrf = CSRFProtect()
app = Flask(__name__, static_folder="./static", template_folder="./templates")
app.config["SECRET_KEY"] = "#nksdhYVi8yGY^r5ryUIHBv"
HOST = "localhost"
PORT = 8000
csrf.init_app(app)

@app.route('/home')
def home():
    return render_template("home.html")

@app.route('/processExperimentalData', methods=['POST'])
def processExperimentalData():
    print(request.form)
    return jsonify({"status": "received"})

if __name__ == '__main__':
    app.run(host=HOST, port=PORT, debug=True)