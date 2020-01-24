# API Definition

Version 0.1

Overview of the input/output requirements for VECTOR. Initially, there will be a single POST endpoint for the frontend to submit parameters to. The backend will then parse the payload and run VECTOR with the proper input files to get a return. The return from the VECTOR code will then be put in the response payload to send back to the frontend.

## Inputs

objectType: string

    sphere, cylinder, plate, plate model

diameter: float

    The diameter of the object [m]

length: float

    The length of the object [m]

area: float

    The area of the object [m^2]

pitch: float

    Pitch angle of the object [deg]

sideslip: float

    Sideslip angle of the object [deg]

temperature: float

    Ambient temperature of the atmosphere [K]

speed: float

    Speed of the object [m/s]

composition: dictionary

    Number density composition of atmospheric consituents [/m3]

    O: float

    O2: float

    N2: float

    He: float

    H: float

accomodationModel: string

    SESAM, Goodman, Fixed

accomodationParameters: dictionary

    Parameters to pass into the accomodation model

    alpha: float

    ms: float

        [amu]

### Example Input Payload

{
"objectType": "cylinder",
"diameter": 1.25,
"length": 2.5,
"area": 0,
"pitch": 30,
"sideslip": 0,
"temperature": 800,
"speed": 7750,
"composition": {"O": 1e11,
                "O2": 1e6,
                "N2": 1e6,
                "He": 1e6,
                "H": 1e4},
"accomodationModel": "SESAM",
"accomodationParameters": {"alpha": 0,
                            "ms": 0}
}

## Outputs

dragCoefficient: float

    The coefficient of drag of the object [-]

projectedArea: float

    The projected area of the object

forceCoefficent: float

    The force coefficient of the object

### Example Output Payload

{
"dragCoefficient": 2.8,
"projectedArea": 3,
"forceCoefficient": 1.2
}
