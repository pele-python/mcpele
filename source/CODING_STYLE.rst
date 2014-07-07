For this project, as well as pele, we use the K&R c++ coding standard (more or
less) as a default

http://en.wikipedia.org/wiki/Indent_style

we use 4 spaces for indentation and *no tabs*

In generall, this means it looks like this

    int main(int argc, char *argv[])
    {
        ...
        while (x == y) {
            something();
            somethingelse();
     
            if (some_error) {
                do_correct();
            } else {
                continue_as_usual();
            }
        }
     
        finalthing();
        ...
    }

For a class it would be something like this

    class MyClass{
    protected:
        double v1;
        double v2;
    public:
        MyClass()
            : v1(0),
              v2(0)
        { 
            do_something();
        }
    };

Some big things I noticed that we need to fix are

if and while statements have the brace on the same line

    if (blah) {
        do something
    }

functions have braces on separate lines

    void func()
    {
        do something
    }

always add white space after commas and around most operators (this is a pet peeve of mine ;)) )

    func(var1,var2,var3);   // no!
    func(var1, var2, var3)  // yes!
    a=b+4; //no!
    a = b + 4; //yes!

with initializer lists, put the colon on a separate line from the name.  And the braces also

    Minimizer::Minimizer()
        : val(0),
          ptr(1)
    {
        do something
    }

Try to keep the lines not too much longer than 80 characters.

Try to generally use braces with if statements.
for loops should never not have braces.  it's just too dangerous.

    // very easy to introduce problems
    if (condition)
        do_something;


add a space after for, if, while, etc

    for(i = 0; i < N; ++i){ //no!
    for (i = 0; i < N; ++i) { //yes!

Put whitespace between operators

    // way too little whitespace
    std::vector<energy_t> energies()const{return property_listing<energy_t>(&Minimum::energy);}
    // better, but it's still quite long and hard to read
    std::vector<energy_t> energies() const { return property_listing<energy_t>(&Minimum::energy); }
    // best
    std::vector<energy_t> energies() const 
    { 
        return property_listing<energy_t>(&Minimum::energy); 
    }


