
int outputPin = 9;
String incomingByte;
int inputNumber;

void setup() {
  Serial.begin(9600);
  pinMode(outputPin, OUTPUT);
}

void loop() {

  if (Serial.available() > 0) {
    
      incomingByte = Serial.readString();
      inputNumber = incomingByte.toInt();
      
      Serial.print("I received: ");      
      Serial.println(inputNumber);
      
      if((inputNumber >= 0) && (inputNumber <= 255)){
        analogWrite(outputPin, inputNumber);   
        Serial.println("Sending");
      }else{
        Serial.println("Error: The number should be 0<number<255");    
      }
  }
}

