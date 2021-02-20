package edu.washington.gs.maccoss.encyclopedia.utils.massspec;

public abstract class Polymer {
	
	public abstract String getName();
	public abstract double getMass(int n);
	public abstract byte getCharge();

	public static Polymer peg=new Polymer() {
		public String getName() {
			return "PEG";
		}
		public double getMass(int n) {
			//[C2H4O]nH2O
			return MassConstants.oh2+n*(MassConstants.carbonMass*2+MassConstants.hydrogenMass*4+MassConstants.oxygenMass)+MassConstants.protonMass;
		}
		
		@Override
		public byte getCharge() {
			return 1;
		}
	};

	public static Polymer pegP2=new Polymer() {
		public String getName() {
			return "PEG+2H";
		}
		public double getMass(int n) {
			//[C2H4O]nH2O
			return (MassConstants.oh2+n*(MassConstants.carbonMass*2+MassConstants.hydrogenMass*4+MassConstants.oxygenMass)+2.0*MassConstants.protonMass)/2.0;
		}
		
		@Override
		public byte getCharge() {
			return 2;
		}
	};

	public static Polymer pegP3=new Polymer() {
		public String getName() {
			return "PEG+3H";
		}
		public double getMass(int n) {
			//[C2H4O]nH2O
			return (MassConstants.oh2+n*(MassConstants.carbonMass*2+MassConstants.hydrogenMass*4+MassConstants.oxygenMass)+3.0*MassConstants.protonMass)/3.0;
		}
		
		@Override
		public byte getCharge() {
			return 3;
		}
	};
	
	public static Polymer ppg=new Polymer() {
		public String getName() {
			return "PPG";
		}
		public double getMass(int n) {
			//[C3H6O]nH2O
			return MassConstants.oh2+n*(MassConstants.carbonMass*3+MassConstants.hydrogenMass*6+MassConstants.oxygenMass)+MassConstants.protonMass;
		}
		
		@Override
		public byte getCharge() {
			return 1;
		}
	};
	
	public static Polymer triton=new Polymer() {
		public String getName() {
			return "Triton";
		}
		public double getMass(int n) {
			// C14H22O(C2H4O)n
			return MassConstants.carbonMass*14+MassConstants.hydrogenMass*22+MassConstants.oxygenMass
					+n*(MassConstants.carbonMass*2+MassConstants.hydrogenMass*4+MassConstants.oxygenMass)+MassConstants.protonMass;
		}
		
		@Override
		public byte getCharge() {
			return 1;
		}
	};
	
	public static Polymer reducedTriton=new Polymer() {
		public String getName() {
			return "TritonR";
		}
		public double getMass(int n) {
			// C14H28O(C2H4O)n 
			return MassConstants.carbonMass*14+MassConstants.hydrogenMass*28+MassConstants.oxygenMass
					+n*(MassConstants.carbonMass*2+MassConstants.hydrogenMass*4+MassConstants.oxygenMass)+MassConstants.protonMass;
		}
		
		@Override
		public byte getCharge() {
			return 1;
		}
	};
	
	public static Polymer tritonNa=new Polymer() {
		public String getName() {
			return "Triton-Na";
		}
		public double getMass(int n) {
			// C14H22O(C2H4O)n
			return MassConstants.carbonMass*14+MassConstants.hydrogenMass*22+MassConstants.oxygenMass
					+n*(MassConstants.carbonMass*2+MassConstants.hydrogenMass*4+MassConstants.oxygenMass)+MassConstants.sodiumMass;
		}
		
		@Override
		public byte getCharge() {
			return 1;
		}
	};
	
	public static Polymer reducedTritonNa=new Polymer() {
		public String getName() {
			return "TritonR-Na";
		}
		public double getMass(int n) {
			// C14H28O(C2H4O)n 
			return MassConstants.carbonMass*14+MassConstants.hydrogenMass*28+MassConstants.oxygenMass
					+n*(MassConstants.carbonMass*2+MassConstants.hydrogenMass*4+MassConstants.oxygenMass)+MassConstants.sodiumMass;
		}
		
		@Override
		public byte getCharge() {
			return 1;
		}
	};
	
	public static Polymer tritonX101=new Polymer() {
		public String getName() {
			return "TritonX-101";
		}
		public double getMass(int n) {
			// C15H24O(C2H4O)n  
			return MassConstants.carbonMass*15+MassConstants.hydrogenMass*24+MassConstants.oxygenMass
					+n*(MassConstants.carbonMass*2+MassConstants.hydrogenMass*4+MassConstants.oxygenMass)+MassConstants.protonMass;
		}
		
		@Override
		public byte getCharge() {
			return 1;
		}
	};
	
	public static Polymer tritonX101R=new Polymer() {
		public String getName() {
			return "TritonX-101R";
		}
		public double getMass(int n) {
			// C15H30O(C2H4O)n 
			return MassConstants.carbonMass*15+MassConstants.hydrogenMass*30+MassConstants.oxygenMass
					+n*(MassConstants.carbonMass*2+MassConstants.hydrogenMass*4+MassConstants.oxygenMass)+MassConstants.protonMass;
		}
		
		@Override
		public byte getCharge() {
			return 1;
		}
	};
	
	public static Polymer polysiloxane=new Polymer() {
		public String getName() {
			return "Polysiloxane";
		}
		public double getMass(int n) {
			//  (C2H6SiO)n
			return n*(MassConstants.carbonMass*2+MassConstants.hydrogenMass*6+MassConstants.siliconMass+MassConstants.oxygenMass)+MassConstants.protonMass;
		}
		
		@Override
		public byte getCharge() {
			return 1;
		}
	};
	
	public static Polymer tween20=new Polymer() {
		public String getName() {
			return "Tween-20";
		}
		public double getMass(int n) {
			// [C18H34O6][C2H4O]n
			return MassConstants.carbonMass*18+MassConstants.hydrogenMass*34+MassConstants.oxygenMass*6
					+n*(MassConstants.carbonMass*2+MassConstants.hydrogenMass*4+MassConstants.oxygenMass)+MassConstants.sodiumMass;
		}
		
		@Override
		public byte getCharge() {
			return 1;
		}
	};
	
	public static Polymer tween40=new Polymer() {
		public String getName() {
			return "Tween-40";
		}
		public double getMass(int n) {
			// [C22H42O6][C2H4O]n 
			return MassConstants.carbonMass*22+MassConstants.hydrogenMass*42+MassConstants.oxygenMass*6
					+n*(MassConstants.carbonMass*2+MassConstants.hydrogenMass*4+MassConstants.oxygenMass)+MassConstants.sodiumMass;
		}
		
		@Override
		public byte getCharge() {
			return 1;
		}
	};
	
	public static Polymer tween60=new Polymer() {
		public String getName() {
			return "Tween-60";
		}
		public double getMass(int n) {
			// [C24H46O6][C2H4O]n
			return MassConstants.carbonMass*24+MassConstants.hydrogenMass*46+MassConstants.oxygenMass*6
					+n*(MassConstants.carbonMass*2+MassConstants.hydrogenMass*4+MassConstants.oxygenMass)+MassConstants.sodiumMass;
		}
		
		@Override
		public byte getCharge() {
			return 1;
		}
	};
	
	public static Polymer tween80=new Polymer() {
		public String getName() {
			return "Tween-80";
		}
		public double getMass(int n) {
			// [C24H44O6][C2H4O]n
			return MassConstants.carbonMass*24+MassConstants.hydrogenMass*44+MassConstants.oxygenMass*6
					+n*(MassConstants.carbonMass*2+MassConstants.hydrogenMass*4+MassConstants.oxygenMass)+MassConstants.sodiumMass;
		}
		
		@Override
		public byte getCharge() {
			return 1;
		}
	};
	public static Polymer[] polymers=new Polymer[] {peg, pegP2, pegP3, ppg, triton, reducedTriton, tritonNa, reducedTritonNa, tritonX101, tritonX101R, polysiloxane, tween20, tween40, tween60, tween80};
}
