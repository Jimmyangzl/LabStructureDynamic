class basic_math:
    
    def __init__(self,a,b=1):
        self.a = a 
        self.b = b
        
    def add(self):
        return self.a+self.b
    
    def sub(self):
        return self.a-self.b
    
    def mul(self):
        return self.a*self.b
    
    def div(self):
        return self.a/self.b
    
    def res(self):
        print('Add: ' + str(self.add())+ '. Sub: '+str(self.sub()) + 
              '. Mul: '+str(self.mul()) + '. Div: '+str(self.div()) +'.' )    