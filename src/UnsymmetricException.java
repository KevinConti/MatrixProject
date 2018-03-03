public class UnsymmetricException extends Exception{
    // Parameterless Constructor
    public UnsymmetricException(
    ) {super("The matrix cannot use Jacobi's, as it is not symmetric");}

    // Constructor that accepts a message
    public UnsymmetricException(String message)
    {
        super(message);
    }
}
