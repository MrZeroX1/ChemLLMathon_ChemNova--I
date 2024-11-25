import os
from flask import Flask, render_template, request, send_from_directory
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image


app = Flask(__name__)

# Configure the upload folder
UPLOAD_FOLDER = os.path.join('static', 'images')
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

# Ensure the upload folder exists
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

# Route for the home page
@app.route('/')
def index():
    return render_template('index.html')

# Route for the about page
@app.route('/about')
def about():
    return render_template('about.html')

# Route for the contact page
@app.route('/contact')
def contact():
    return render_template('contact.html')

@app.route('/generate', methods=['POST'])
def generate():
    # Get JSON data from the request
    data = request.get_json()
    smiles = data.get('smiles')  # Extract the SMILES string

    if smiles:
        try:
            # Create a molecule object from the SMILES string
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # Generate a unique filename for the molecule image
                filename = f"molecule_{hash(smiles)}.png"
                image_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)

                # Use RDKit's Draw method with transparent background
                d2d = rdMolDraw2D.MolDraw2DCairo(300, 300)  # 300x300 image
                d2d.drawOptions().bgColor = None  # Transparent background
                d2d.DrawMolecule(mol)
                d2d.FinishDrawing()

                # Save the generated molecule image
                with open(image_path, "wb") as img_file:
                    img_file.write(d2d.GetDrawingText())

                # Return the image path as a JSON response
                return {'image_path': f"images/{filename}"}, 200
            else:
                # Handle invalid SMILES string
                return {'error': 'Invalid SMILES string. Please try again.'}, 400
        except Exception as e:
            # Handle errors during image generation
            return {'error': f"An error occurred: {str(e)}"}, 500
    else:
        # Handle empty input
        return {'error': 'No SMILES string provided. Please enter a valid SMILES.'}, 400
    
if __name__ == '__main__':
    app.run(debug=True)
