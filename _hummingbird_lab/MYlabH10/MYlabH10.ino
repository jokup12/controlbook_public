#include "arduino_setup.h"
#include "time_utilities.h"
#include "sensor_utilities.h"
#include "motor_utilities.h"
#include "reference_utilities.h"
#include "ctrlPID.h"

//=============================================================================
// declare global structures
//=============================================================================

// instantiate structures
Reference reference;

// instantiate classes as global variables
TimeUtilities timing;
SensorUtilities sensors;
MotorUtilities rotors;
ReferenceUtilities signal_generator;
CtrlPID controller;

//=============================================================================
// arduino setup function (runs once at start of simulation)
//=============================================================================
void setup()
{
  // start serial communication
  Serial.begin(9600);

  // initialize all classes
  timing.init();  // initialize current time and sample rate
  sensors.init();  // initialize sensors
  controller.init();  // initialize controller
  

  // initialize buttons on breakout board and watchdog timers
  initialize_buttons();
}

//=============================================================================
// arduino loop function (loops forever)
//=============================================================================
void loop()
{
  timing.update();  // update current time and sample rate
  sensors.update();  // update sensors
    //float psi_ref = signal_generator.joystick_yaw();

  signal_generator.init(15*3.14/180, 0.05, 0.0);
  float psi_ref = signal_generator.square_signal(timing.current);
  //signal_generator.init(15*3.14/180, 0.1, 0.0); 
  float theta_ref = signal_generator.square_signal(timing.current);

  controller.update(psi_ref, theta_ref, sensors, rotors, timing.Ts);  // update controller

  // zero encoders if zero button pushed
  zeroButton.update();  
  if ( zeroButton.pressed() ) sensors.zero(); 
   
  // calibrate rotors if calibrate button pushed
  calButton.update();
  if ( calButton.pressed() ) rotors.init();
  
  // reset watchdog timer
  wdt_reset();
  digitalWrite(LED_RX, LOW);
}
