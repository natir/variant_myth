//! A variant annotater.

#![warn(missing_docs)]

/* std use */

/* crate use */

/* project use */


/* mod declaration */
pub mod cli;
pub mod error;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
