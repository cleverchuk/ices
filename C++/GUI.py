## GUI for finding maximum feasible combinations

from Tkinter import *
import subprocess



class GUI(Tk):
   
    
    def __init__(self, binary):
        Tk.__init__(self)
        self.bin = "./"+binary
        # value holders
        self.A = IntVar()
        self.B = IntVar()
        self.GD = IntVar()
        self.GE = IntVar()
        self.GM = IntVar()
        self.SE = IntVar()
        self.SD = IntVar()
        self.VD = IntVar()

        self.time = StringVar()
        self.params =[self.A,self.B,self.GD,self.GE,self.GM,
                      self.SD,self.SE,self.VD]

        
        self.container = Frame(self)
        self.container.pack(fill=BOTH, expand=1)

        # widgets
        cb_A = Checkbutton(self.container,text="A", variable = self.A,
                              onvalue="1")
        cb_B = Checkbutton(self.container, variable = self.B,text="B",
                              onvalue="2")

        cb_GD = Checkbutton(self.container, variable = self.GD,text="Grainboundary diffusion",
                              onvalue="3")
        cb_GE = Checkbutton(self.container, variable = self.GE,text="Grainboundary energy",
                              onvalue="4")
        
        cb_GM = Checkbutton(self.container, variable = self.GM,text="Grainboundary mobility",
                              onvalue="5")
        cb_SD = Checkbutton(self.container, variable = self.SD,text="Surface diffusion",
                              onvalue="6")
        
        cb_SE = Checkbutton(self.container, variable = self.SE,text="Surface energy",
                              onvalue="7")
        cb_VD = Checkbutton(self.container, variable = self.VD,text="Volume diffusion",
                              onvalue="8")

        simulTime = Entry(self.container, width=5,text="final timer", textvariable=self.time)
        simulTime.grid(row=9, column=2, sticky="WE")

        label = Label(self.container,text="Final time",bg="black",fg="white")
        label.grid(row=9,column=1,sticky="W")
                
        command = Button(self.container,bg="green",text="Run",
                         command=self.process, width=10)
        command.grid(row=10,column=2, sticky="E")
        
        

        cb_A.grid(row=0,column=0, sticky="W")
        cb_B.grid(row=1,column=0, sticky="W")
        cb_GD.grid(row=3,column=0, sticky="W")
        cb_GE.grid(row=2,column=0, sticky="W")
        cb_GM.grid(row=4,column=0, sticky="W")
        cb_SD.grid(row=5,column=0, sticky="W")
        cb_SE.grid(row=6,column=0, sticky="W")
        cb_VD.grid(row=7,column=0, sticky="W")


    def process(self):
        temp = list()
        cmd = []
        for var in self.params:
            value = var.get()
            if 0 != value:
                value -= 1
                temp.append(str(value))

        time = self.time.get()
        num = len(temp)
        cmd.append(self.bin)
        
        cmd.append(time)
        cmd.append(time)
        cmd.append(str(num))

        
        cmd = cmd + temp
        subprocess.call(cmd)
                    
                    

if __name__ == "__main__":
    app = GUI("combinations") # Enter the name of the binary in the parenthesis
    app.geometry("350x250")
    app.title("Parameter search App")
    app.mainloop()
