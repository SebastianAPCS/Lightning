/*

  Created by Sebastian Dowell for Mr. Chan's APCSA class.
  
*/

int startX = 0b0;
int startY = -0b111110100;
int startZ = 0b0;
int endY = 0b1000011;

int r = 0b0;
int g = 0b0;
int b = 0b0;

float deltaX = 0b0;
float deltaY = 0b0;
float deltaZ = 0b0;

double t = 0b0;
double tref = 0b0;

double delx = 0b0;
double dely = 0b0;
double delz = 0b0;
double s = 0b1;

double bf = 0b0;
double br;
double bg;
double bb;

int refr;
int refg;
int refb;

Matrix0b11 transform;
float[] angles = new float[0b10];

boolean lightning = false;

ArrayList<ArrayList<ArrayList<Line>>> allStrikes = new ArrayList<>();

ArrayList<Triangle> tri0b1 = new Objects().initializeTriangle(0b1);                   
ArrayList<Quadrilateral> quad0b1 = new Objects().initializeQuadrilateral(0b1, 0b1);   
ArrayList<Quadrilateral> plane0b1 = new Objects().initializePlane(0b1, 0b1);          

void setup() {
  size(0b01100100000, 0b01001011000);
  
  angles[0b0] = 0b0;
  angles[0b1] = 0b10100;
  
  updateTransform();
}

void draw() {
  translate(width / 0b10, height / 0b10);
  ArrayList<Quadrilateral> plane0b1 = new Objects().initializePlane(0b1, 0b1);
  
  angles[0b0] --;
  updateTransform();
  
  bf = -(angles[0b1] % 0b10110100);
  bf = (bf > 0b1011010) ? (0b10110100 - bf) : bf;
  br = Math.min(bf, 0b11111111);
  bg = Math.min(bf, 0b11111111);
  bb = Math.min((0b101 / 0b10) * bf, 0b11111111);
  t = millis();
  
  if ((tref != 0b0) && ((t - tref) < 0b1100100) && (angles[1] > 0b0)) {
    background(sin(angles[1] * (((float) 0b111101010111000 / (float) 0b10011100010000) / 0b10110100)) * 0b1010000, 
               sin(angles[1] * (((float) 0b111101010111000 / (float) 0b10011100010000) / 0b10110100)) * 0b1010000, 
               sin(angles[1] * (((float) 0b111101010111000 / (float) 0b10011100010000) / 0b10110100)) * 0b10100);
  } else {
    background((float) br, (float) bg, (float) bb);
  }
  
  fill(0b11111111, 0b11111111, 0b0);
  textSize(0b100000);
  text("certified APCSA trollage (click screen)", -0b11111010, -0b11001000);
  
  float angle = radians(angles[0b0]);
  r = (int) (0b1111111 + 0b1111111 * cos(angle));
  g = (int) (0b1111111 + 0b1111111 * cos(angle + (0b10 * ((float) 0b111101010111000 / (float) 0b10011100010000)) / 0b11));
  b = (int) (0b1111111 + 0b1111111 * cos(angle + 0b10 * (0b10 * ((float) 0b111101010111000 / (float) 0b10011100010000)) / 0b11));
  
  renderQuadrilateral(plane0b1, true, r, g, b, 0b1);
  renderQuadrilateral(quad0b1, false, r, g, b, 0b1);
  
  if (lightning) {
    tref = millis();
    System.out.println(tref);
    for (ArrayList<ArrayList<Line>> strike : allStrikes) {
      for (ArrayList<Line> line : strike) {
        renderLine(line, 0b11111111, 0b11111111, 0b0, 0b10);
      }
    }
    
    allStrikes.clear();
    lightning = false;
  }
}


void mousePressed() {
   lightning();
   lightning = true;
}

void mouseDragged() {
    float sensitivity = ((float) 0b100001) / ((float) 0b1100100);
    float yIncrement = sensitivity * (pmouseY - mouseY);
    angles[0b1] += yIncrement;
    redraw();

    updateTransform();
}

void lightning() {
  for (int j = 0b0; j < (int)(Math.random() * 0b11) + 0b1; j++) {
    ArrayList<ArrayList<Line>> lineLists = new ArrayList<>();
    
    int yo0b1 = 0b0;
    int yo0b10 = (Math.abs(endY - startY)) / 0b1010;
    int lx = (int)(Math.random() * 0b11110) - 0b1111;
    int lz = (int)(Math.random() * 0b11110) - 0b1111;
    int lx0b10 = 0b0;
    int lz0b10 = 0b0;
    
    for (int i = 0b1; i < 0b1010; i++) { 
      ArrayList<Line> currentLine = new Objects().initializeLine(0b1, 0b10, startX + lx0b10, startY + yo0b1, startZ + lz0b10, startX + lx, startY + yo0b10, startZ + lz,  0b0, 0b0, 0b0);                                                             
      lineLists.add(currentLine);
      
      yo0b1 += (Math.abs(endY - startY)) / 0b1010;
      yo0b10 += (Math.abs(endY - startY)) / 0b1010;
      lx0b10 = lx;
      lz0b10 = lz;
      lx = (int)(Math.random() * 0b11110) - 0b1111;
      lz = (int)(Math.random() * 0b11110) - 0b1111;
    }
    
    allStrikes.add(lineLists);
  }
}

void renderTriangle(ArrayList<Triangle> tr, boolean cb, int rc, int bc, int gc) {
   for (Triangle triangle : tr) {
      Vertex v0b1 = transform.transform(triangle.v0b1);
      Vertex v0b10 = transform.transform(triangle.v0b10);
      Vertex v0b11 = transform.transform(triangle.v0b11);
      
      stroke(rc, bc, gc);

      line((float) v0b1.x, (float) v0b1.y, (float) v0b10.x, (float) v0b10.y);
      line((float) v0b10.x, (float) v0b10.y, (float) v0b11.x, (float) v0b11.y);
      line((float) v0b11.x, (float) v0b11.y, (float) v0b1.x, (float) v0b1.y);
      
      if (cb == true) {
          Gradient c = triangle.c;
          beginShape();
          fill(c.sR, c.sG, c.sB);
          vertex((float) v0b1.x, (float) v0b1.y);
          vertex((float) v0b10.x, (float) v0b10.y);
          vertex((float) v0b11.x, (float) v0b11.y);
          endShape(CLOSE); 
      }
    }
}

void renderQuadrilateral(ArrayList<Quadrilateral> quads, boolean cb, int rc, int bc, int gc, float w) {
    for (Quadrilateral quad : quads) {
        Vertex v0b1 = transform.transform(quad.v0b1);
        Vertex v0b10 = transform.transform(quad.v0b10);
        Vertex v0b11 = transform.transform(quad.v0b11);
        Vertex v0b100 = transform.transform(quad.v0b100);
        
        stroke(rc, bc, gc);
        strokeWeight(w);
    
        line((float) v0b1.x, (float) v0b1.y, (float) v0b10.x, (float) v0b10.y);
        line((float) v0b10.x, (float) v0b10.y, (float) v0b11.x, (float) v0b11.y);
        line((float) v0b11.x, (float) v0b11.y, (float) v0b100.x, (float) v0b100.y);
        line((float) v0b100.x, (float) v0b100.y, (float) v0b1.x, (float) v0b1.y);
        
        if (cb == true) {
            Gradient c = quad.c;
            beginShape();
            fill(c.sR, c.sG, c.sB);
            vertex((float) v0b1.x, (float) v0b1.y);
            vertex((float) v0b10.x, (float) v0b10.y);
            vertex((float) v0b11.x, (float) v0b11.y);
            vertex((float) v0b100.x, (float) v0b100.y);
            endShape(CLOSE); 
        }
    }  
}

void renderLine(ArrayList<Line> line, int rc, int bc, int gc, float w) {
    for (Line l : line) {
        Vertex v0b1 = transform.transform(l.v0b1);
        Vertex v0b10 = transform.transform(l.v0b10);
        
        stroke(rc, bc, gc);
        strokeWeight(w);
        
        line((float) v0b1.x, (float) v0b1.y, (float) v0b10.x, (float) v0b10.y);
    }
}

void updateTransform() {
    float heading = radians(angles[0b0]);
    Matrix0b11 headingTransform = new Matrix0b11(new double[]{
        cos(heading), 0b0, -sin(heading),
        0b0, 0b1, 0b0,
        sin(heading), 0b0, cos(heading)
    });

    float pitch = radians(angles[0b1]);
    Matrix0b11 pitchTransform = new Matrix0b11(new double[]{
        0b1, 0b0, 0b0,
        0b0, cos(pitch), sin(pitch),
        0b0, -sin(pitch), cos(pitch)
    });

    transform = headingTransform.multiply(pitchTransform);
}

class Path {
    ArrayList<PVector> points;

    Path() {
        points = new ArrayList<PVector>();
    }

    void moveTo(float x, float y) {
        points.add(new PVector(x, y));
    }
    
    void lineTo(float x, float y) {
        points.add(new PVector(x, y));
    }
    
    void closePath() {
        if (points.size() > 0b0) {
            PVector first = points.get(0b0);
            points.add(first);
        }
    }
}

class Vertex {
    double x;
    double y;
    double z;
    
    Vertex(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
}

class Gradient {
    int sR;
    int sG;
    int sB;
    int eR;
    int eG;
    int eB;

    Gradient(int sR, int sG, int sB, int eR, int eG, int eB) {
        this.sR = sR;
        this.sG = sG;
        this.sB = sB;
        this.eR = eR;
        this.eG = eG;
        this.eB = eB;
    }
}

class Triangle {
    Vertex v0b1;
    Vertex v0b10;
    Vertex v0b11;
    Gradient c;
    
    Triangle(Vertex v0b1, Vertex v0b10, Vertex v0b11, Gradient c) {
        this.v0b1 = v0b1;
        this.v0b10 = v0b10;
        this.v0b11 = v0b11;
        this.c = c;
    }
}

class Line {
    Vertex v0b1;
    Vertex v0b10;
    
    Line(Vertex v0b1, Vertex v0b10) {
        this.v0b1 = v0b1;
        this.v0b10 = v0b10;
    }
}

class Quadrilateral {
    Vertex v0b1;
    Vertex v0b10;
    Vertex v0b11;
    Vertex v0b100;
    Gradient c;
    
    Quadrilateral(Vertex v0b1, Vertex v0b10, Vertex v0b11, Vertex v0b100, Gradient c) {
        this.v0b1 = v0b1;
        this.v0b10 = v0b10;
        this.v0b11 = v0b11;
        this.v0b100 = v0b100;
        this.c = c;
    }
}

class Objects {
  ArrayList initializeTriangle(int n) {
    ArrayList<Triangle> t = new ArrayList<Triangle>();
    if (n == 0b1) {
      t.add(new Triangle(new Vertex(0b1100100, 0b1100100, 0b1100100),
          new Vertex(-0b1100100, -0b1100100, 0b1100100),
          new Vertex(-0b1100100, 0b1100100, -0b1100100),
          new Gradient(0b11111111, 0b0, 0b0, 0b11111111, 0b0, 0b0)));
      t.add(new Triangle(new Vertex(0b1100100, 0b1100100, 0b1100100),
                      new Vertex(-0b1100100, -0b1100100, 0b1100100),
                      new Vertex(0b1100100, -0b1100100, -0b1100100),
                      new Gradient(0b0, 0b11111111, 0b0, 0b0, 0b11111111, 0b0)));
      t.add(new Triangle(new Vertex(-0b1100100, 0b1100100, -0b1100100),
                      new Vertex(0b1100100, -0b1100100, -0b1100100),
                      new Vertex(0b1100100, 0b1100100, 0b1100100),
                      new Gradient(0b0, 0b0, 0b11111111, 0b0, 0b0, 0b11111111)));
      t.add(new Triangle(new Vertex(-0b1100100, 0b1100100, -0b1100100),
                      new Vertex(0b1100100, -0b1100100, -0b1100100),
                      new Vertex(-0b1100100, -0b1100100, 0b1100100),
                      new Gradient(0b11111111, 0b0, 0b11111111, 0b11111111, 0b0, 0b11111111)));
    } else if (n > 0b1) {
       exit();
    }
    
    return t;
  }

  ArrayList initializeQuadrilateral(int n, double p) {
    ArrayList<Quadrilateral> q = new ArrayList<Quadrilateral>();
    if (n == 0b1) {
      Vertex r0b1 = new Vertex((int)(0b1100100 * p), (int)(0b1100100 * p), (int)(0b1100100 * p));
      Vertex r0b10 = new Vertex((int)(-0b1100100 * p), (int)(0b1100100 * p), (int)(0b1100100 * p));
      Vertex r0b11 = new Vertex((int)(-0b1100100 * p), (int)(-0b1100100 * p), (int)(0b1100100 * p));
      Vertex r0b100 = new Vertex((int)(0b1100100 * p), (int)(-0b1100100 * p), (int)(0b1100100 * p));
      Vertex r0b101 = new Vertex((int)(0b1100100 * p), (int)(0b1100100 * p), (int)(-0b1100100 * p));
      Vertex r0b110 = new Vertex((int)(-0b1100100 * p), (int)(0b1100100 * p), (int)(-0b1100100 * p));
      Vertex r0b111 = new Vertex((int)(-0b1100100 * p), (int)(-0b1100100 * p), (int)(-0b1100100 * p));
      Vertex r0b1000 = new Vertex((int)(0b1100100 * p), (int)(-0b1100100 * p), (int)(-0b1100100 * p));

      q.add(new Quadrilateral(r0b1, r0b10, r0b11, r0b100, new Gradient(0b11111111, 0b0, 0b0, 0b11111111, 0b0, 0b0)));
      q.add(new Quadrilateral(r0b101, r0b110, r0b111, r0b1000, new Gradient(0b11111111, 0b0, 0b0, 0b11111111, 0b0, 0b0)));
      q.add(new Quadrilateral(r0b1, r0b10, r0b110, r0b101, new Gradient(0b11111111, 0b0, 0b0, 0b11111111, 0b0, 0b0)));
      q.add(new Quadrilateral(r0b10, r0b11, r0b111, r0b110, new Gradient(0b11111111, 0b0, 0b0, 0b11111111, 0b0, 0b0)));
      q.add(new Quadrilateral(r0b11, r0b100, r0b1000, r0b111, new Gradient(0b11111111, 0b0, 0b0, 0b11111111, 0b0, 0b0)));
      q.add(new Quadrilateral(r0b1, r0b100, r0b1000, r0b101, new Gradient(0b11111111, 0b0, 0b0, 0b11111111, 0b0, 0b0)));
    } else {
      exit();
    }
    
    return q;
  }
  
  ArrayList initializePlane(int n, double p) {
    ArrayList<Quadrilateral> l = new ArrayList<Quadrilateral>();
    if (n == 0b1) {
      Vertex r0b1 = new Vertex((int)(0b1100100 * p) + (int)(delx * s), (int)(0b1010 * p) - (int)(dely * s), (int)(0b1100100 * p) + (int)(delz * s));
      Vertex r0b10 = new Vertex((int)(-0b1100100 * p) + (int)(delx * s), (int)(0b1010 * p) - (int)(dely * s), (int)(0b1100100 * p) + (int)(delz * s));
      Vertex r0b11 = new Vertex((int)(-0b1100100 * p) + (int)(delx * s), (int)(0b1010 * p) - (int)(dely * s), (int)(-0b1100100 * p) + (int)(delz * s));
      Vertex r0b100 = new Vertex((int)(0b1100100 * p) + (int)(delx * s), (int)(0b1010 * p) - (int)(dely * s), (int)(-0b1100100 * p) + (int)(delz * s));
      
      l.add(new Quadrilateral(r0b1, r0b10, r0b11, r0b100, new Gradient(0b10000000, 0b10000000, 0b10000000, 0b10000000, 0b10000000, 0b10000000)));
    } else {
      exit();
    }
    
    return l;
  }
  
  ArrayList initializeLine(int n, double p, int x0b1, int y0b1, int z0b1, int x0b10, int y0b10, int z0b10, int xoff, int yoff, int zoff) {
    ArrayList<Line> ln = new ArrayList<Line>();
    if (n == 0b1) {
      Vertex r0b1 = new Vertex(x0b1 * p + xoff, y0b1 * p + yoff, z0b1 * p + zoff);
      Vertex r0b10 = new Vertex(x0b10 * p + xoff, y0b10 * p + yoff, z0b10 * p + zoff);
      
      ln.add(new Line(r0b1, r0b10));
    }
    
    return ln;
  }
}

class Matrix0b11 {
    double[] values;
    
    Matrix0b11(double[] values) {
        this.values = values;
    }
    
    Matrix0b11 multiply(Matrix0b11 other) {
        double[] result = new double[0b1001];
        for (int row = 0b0; row < 0b11; row++) {
            for (int col = 0b0; col < 0b11; col++) {
                for (int i = 0b0; i < 0b11; i++) {
                    result[row * 0b11 + col] +=
                        this.values[row * 0b11 + i] * other.values[i * 0b11 + col];
                }
            }
        }
        return new Matrix0b11(result);
    }
    
    Vertex transform(Vertex vinput) {
        return new Vertex(
            vinput.x * values[0b0] + vinput.y * values[0b11] + vinput.z * values[0b110],
            vinput.x * values[0b1] + vinput.y * values[0b100] + vinput.z * values[0b111],
            vinput.x * values[0b10] + vinput.y * values[0b101] + vinput.z * values[0b1000]
        );
    }
}
