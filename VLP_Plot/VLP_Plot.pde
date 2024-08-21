import processing.serial.*;

Serial teensy;
float triangleRadius = 30;
String line;
int scalingFactor = 10;

void setup()
{
  //size(1920,1080);
  fullScreen();
  String portName = Serial.list()[3]; 
  
  println(portName);
  
  teensy = new Serial(this, portName, 115200);
  frameRate(120);
  
}

void draw()
{
  if ( teensy.available() > 0) {  
    char c = (char) teensy.read();
    if (c != '\n') {
      line += c;
    } else {
      String[] list = split(line, " ");
      float x = float(list[0]);
      float y = float(list[1]);
      float theta = float(list[2]);
      push();
      translate(width/4, height/2);
      rotate(-PI/2);
      //scale(2);
      background(0);
      drawLed(0, 0);
      drawLed(41, 0);
      drawLed(0, 41);
      drawLed(41, 41);
      drawDevice(x*100, y*100, theta);
      pop();

      line = "";
      teensy.clear();
    }
  }
}


void drawLed(float x, float y) {
  rectMode(CENTER);
  fill(255);
  stroke(128);
  rect(x*scalingFactor, y*scalingFactor, 20, 20);
}

void drawDevice(float x, float y, float theta) {
  float x1 = x*scalingFactor + triangleRadius * cos(PI/3 + theta);
  float y1 = y*scalingFactor + triangleRadius * sin(PI/3 + theta);
  float x2 = x*scalingFactor + triangleRadius * cos(PI + theta);
  float y2 = y*scalingFactor + triangleRadius * sin(PI + theta);
  float x3 = x*scalingFactor + triangleRadius * cos(-PI/3 + theta);
  float y3 = y*scalingFactor + triangleRadius * sin(-PI/3 + theta);
  
  fill(255,0,0);
  stroke(128);
  //println(x1,y1);
  triangle(x1, y1, x2, y2, x3, y3);
  stroke(0,255,0);
  line(x1,y1,x3,y3);

}
